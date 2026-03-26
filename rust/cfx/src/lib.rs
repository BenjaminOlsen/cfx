//! Rust bindings for cfx - arbitrary precision integers and cryptographic operations.
//!
//! cfx provides:
//! - **Arbitrary precision integers** (`Big`) with full arithmetic support
//! - **Cryptographic hashes**: SHA-256, BLAKE2b/s
//! - **Digital signatures**: Ed25519
//! - **Key exchange**: X25519
//! - **Random number generation**: ChaCha20-based CSPRNG
//!
//! # Quick Start
//!
//! ```
//! use cfx::{Big, hash::Sha256, ed25519::Ed25519Keypair};
//!
//! // Big integers
//! let a = Big::from(123u64);
//! let b = Big::from(456u64);
//! let c = &a * &b;
//! println!("{} * {} = {}", a, b, c);
//!
//! // Hashing
//! let hash = Sha256::hash(b"hello world");
//!
//! // Signatures
//! let seed = [0u8; 32]; // Use random bytes!
//! let keypair = Ed25519Keypair::from_seed(&seed);
//! let sig = keypair.sign(b"message");
//! assert!(sig.verify(b"message", keypair.public_key()));
//! ```

pub mod ffi;
pub mod error;
pub mod big;
pub mod hash;
pub mod ed25519;
pub mod x25519;
pub mod rand;

// Re-export main types at crate root
pub use big::Big;
pub use error::{Error, Result};
pub use hash::{Sha256, Blake2b, Blake2s};
pub use ed25519::{Ed25519Keypair, Ed25519Signature};
pub use x25519::{X25519PrivateKey, X25519PublicKey, X25519SharedSecret};
pub use rand::ChaChaRng;

#[cfg(test)]
mod tests {
    use super::*;

    // Big integer tests (moved from old lib.rs)

    #[test]
    fn test_new_is_zero() {
        let a = Big::new();
        assert!(a.is_zero());
        assert_eq!(format!("{}", a), "0");
    }

    #[test]
    fn test_from_u64() {
        let a = Big::from(42u64);
        assert_eq!(format!("{}", a), "42");
    }

    #[test]
    fn test_from_u32() {
        let a = Big::from(12345u32);
        assert_eq!(format!("{}", a), "12345");
    }

    #[test]
    fn test_from_usize() {
        let a = Big::from(99999usize);
        assert_eq!(format!("{}", a), "99999");
    }

    #[test]
    fn test_default() {
        let a = Big::default();
        assert!(a.is_zero());
    }

    #[test]
    fn test_clone() {
        let a = Big::from(12345u64);
        let b = a.clone();
        assert_eq!(a, b);
        assert_eq!(format!("{}", b), "12345");
    }

    #[test]
    fn test_debug() {
        let a = Big::from(42u64);
        assert_eq!(format!("{:?}", a), "Big(42)");
    }

    #[test]
    fn test_from_dec() {
        let a = Big::from_dec("123456789012345678901234567890").unwrap();
        assert_eq!(format!("{}", a), "123456789012345678901234567890");
    }

    #[test]
    fn test_from_dec_zero() {
        let a = Big::from_dec("0").unwrap();
        assert!(a.is_zero());
    }

    #[test]
    fn test_from_dec_invalid() {
        assert!(Big::from_dec("not_a_number").is_none());
    }

    #[test]
    fn test_from_hex() {
        let a = Big::from_hex("deadbeef").unwrap();
        assert_eq!(a, Big::from(0xdeadbeefu64));
    }

    #[test]
    fn test_from_hex_large() {
        let a = Big::from_hex("123456789abcdef0123456789abcdef0").unwrap();
        assert_eq!(a.to_hex(), "123456789abcdef0123456789abcdef0");
    }

    #[test]
    fn test_to_hex() {
        let a = Big::from(255u64);
        assert_eq!(a.to_hex(), "ff");
    }

    #[test]
    fn test_to_hex_zero() {
        let a = Big::new();
        assert_eq!(a.to_hex(), "0");
    }

    #[test]
    fn test_is_zero() {
        assert!(Big::new().is_zero());
        assert!(Big::from(0u64).is_zero());
        assert!(!Big::from(1u64).is_zero());
    }

    #[test]
    fn test_is_one() {
        assert!(!Big::new().is_one());
        assert!(Big::from(1u64).is_one());
        assert!(!Big::from(2u64).is_one());
    }

    #[test]
    fn test_is_even_odd() {
        assert!(Big::from(0u64).is_even());
        assert!(Big::from(2u64).is_even());
        assert!(Big::from(100u64).is_even());
        assert!(!Big::from(1u64).is_even());
        assert!(!Big::from(99u64).is_even());

        assert!(Big::from(1u64).is_odd());
        assert!(Big::from(99u64).is_odd());
        assert!(!Big::from(0u64).is_odd());
        assert!(!Big::from(100u64).is_odd());
    }

    #[test]
    fn test_is_prime() {
        assert!(!Big::from(0u64).is_prime());
        assert!(!Big::from(1u64).is_prime());
        assert!(Big::from(2u64).is_prime());
        assert!(Big::from(3u64).is_prime());
        assert!(!Big::from(4u64).is_prime());
        assert!(Big::from(5u64).is_prime());
        assert!(Big::from(7u64).is_prime());
        assert!(!Big::from(9u64).is_prime());
        assert!(Big::from(11u64).is_prime());
        assert!(Big::from(13u64).is_prime());
        assert!(!Big::from(15u64).is_prime());
        assert!(Big::from(17u64).is_prime());
        assert!(Big::from(19u64).is_prime());
        assert!(Big::from(104729u64).is_prime());
        assert!(!Big::from(104730u64).is_prime());
    }

    #[test]
    fn test_bit_len() {
        assert_eq!(Big::new().bit_len(), 0);
        assert_eq!(Big::from(1u64).bit_len(), 1);
        assert_eq!(Big::from(2u64).bit_len(), 2);
        assert_eq!(Big::from(3u64).bit_len(), 2);
        assert_eq!(Big::from(4u64).bit_len(), 3);
        assert_eq!(Big::from(255u64).bit_len(), 8);
        assert_eq!(Big::from(256u64).bit_len(), 9);
    }

    #[test]
    fn test_cmp() {
        let a = Big::from(100u64);
        let b = Big::from(200u64);
        assert!(a < b);
        assert!(b > a);
        assert_eq!(a, a.clone());
    }

    #[test]
    fn test_cmp_eq() {
        let a = Big::from(12345u64);
        let b = Big::from(12345u64);
        assert_eq!(a, b);
        assert!(!(a < b));
        assert!(!(a > b));
    }

    #[test]
    fn test_partial_eq_u64() {
        let a = Big::from(42u64);
        assert_eq!(a, 42u64);
        assert!(a != 43u64);
    }

    #[test]
    fn test_ordering() {
        let a = Big::from(10u64);
        let b = Big::from(20u64);
        let c = Big::from(10u64);
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Less);
        assert_eq!(b.cmp(&a), std::cmp::Ordering::Greater);
        assert_eq!(a.cmp(&c), std::cmp::Ordering::Equal);
    }

    #[test]
    fn test_add() {
        let a = Big::from(100u64);
        let b = Big::from(23u64);
        let c = &a + &b;
        assert_eq!(c, 123u64);
    }

    #[test]
    fn test_add_owned() {
        let a = Big::from(100u64);
        let b = Big::from(23u64);
        let c = a + b;
        assert_eq!(c, 123u64);
    }

    #[test]
    fn test_add_u64() {
        let a = Big::from(100u64);
        let b = a + 50u64;
        assert_eq!(b, 150u64);
    }

    #[test]
    fn test_add_assign() {
        let mut a = Big::from(100u64);
        let b = Big::from(50u64);
        a += &b;
        assert_eq!(a, 150u64);
    }

    #[test]
    fn test_add_assign_u64() {
        let mut a = Big::from(100u64);
        a += 25u64;
        assert_eq!(a, 125u64);
    }

    #[test]
    fn test_sub() {
        let a = Big::from(100u64);
        let b = Big::from(23u64);
        let c = &a - &b;
        assert_eq!(c, 77u64);
    }

    #[test]
    fn test_sub_owned() {
        let a = Big::from(100u64);
        let b = Big::from(30u64);
        let c = a - b;
        assert_eq!(c, 70u64);
    }

    #[test]
    fn test_sub_u64() {
        let a = Big::from(100u64);
        let b = a - 25u64;
        assert_eq!(b, 75u64);
    }

    #[test]
    fn test_sub_assign() {
        let mut a = Big::from(100u64);
        let b = Big::from(40u64);
        a -= &b;
        assert_eq!(a, 60u64);
    }

    #[test]
    fn test_sub_assign_u64() {
        let mut a = Big::from(100u64);
        a -= 35u64;
        assert_eq!(a, 65u64);
    }

    #[test]
    fn test_mul() {
        let a = Big::from(12u64);
        let b = Big::from(34u64);
        let c = &a * &b;
        assert_eq!(c, 408u64);
    }

    #[test]
    fn test_mul_owned() {
        let a = Big::from(7u64);
        let b = Big::from(8u64);
        let c = a * b;
        assert_eq!(c, 56u64);
    }

    #[test]
    fn test_mul_u64() {
        let a = Big::from(25u64);
        let b = a * 4u64;
        assert_eq!(b, 100u64);
    }

    #[test]
    fn test_mul_assign() {
        let mut a = Big::from(10u64);
        let b = Big::from(10u64);
        a *= &b;
        assert_eq!(a, 100u64);
    }

    #[test]
    fn test_mul_assign_u64() {
        let mut a = Big::from(10u64);
        a *= 7u64;
        assert_eq!(a, 70u64);
    }

    #[test]
    fn test_mul_large() {
        let a = Big::from_dec("123456789012345678901234567890").unwrap();
        let b = Big::from(2u64);
        let c = &a * &b;
        assert_eq!(format!("{}", c), "246913578024691357802469135780");
    }

    #[test]
    fn test_div() {
        let a = Big::from(100u64);
        let b = Big::from(7u64);
        let c = &a / &b;
        assert_eq!(c, 14u64);
    }

    #[test]
    fn test_div_owned() {
        let a = Big::from(99u64);
        let b = Big::from(11u64);
        let c = a / b;
        assert_eq!(c, 9u64);
    }

    #[test]
    fn test_div_u64() {
        let a = Big::from(100u64);
        let b = a / 5u64;
        assert_eq!(b, 20u64);
    }

    #[test]
    fn test_div_assign() {
        let mut a = Big::from(100u64);
        let b = Big::from(4u64);
        a /= &b;
        assert_eq!(a, 25u64);
    }

    #[test]
    fn test_div_assign_u64() {
        let mut a = Big::from(100u64);
        a /= 10u64;
        assert_eq!(a, 10u64);
    }

    #[test]
    fn test_rem() {
        let a = Big::from(100u64);
        let b = Big::from(7u64);
        let c = &a % &b;
        assert_eq!(c, 2u64);
    }

    #[test]
    fn test_rem_owned() {
        let a = Big::from(17u64);
        let b = Big::from(5u64);
        let c = a % b;
        assert_eq!(c, 2u64);
    }

    #[test]
    fn test_rem_u64() {
        let a = Big::from(100u64);
        let r = a % 7u64;
        assert_eq!(r, 2u64);
    }

    #[test]
    fn test_rem_u64_ref() {
        let a = Big::from(100u64);
        let r = &a % 7u64;
        assert_eq!(r, 2u64);
    }

    #[test]
    fn test_rem_assign() {
        let mut a = Big::from(100u64);
        let b = Big::from(7u64);
        a %= &b;
        assert_eq!(a, 2u64);
    }

    #[test]
    fn test_pow() {
        let a = Big::from(2u64);
        let b = a.pow(10);
        assert_eq!(b, 1024u64);
    }

    #[test]
    fn test_pow_zero() {
        let a = Big::from(5u64);
        let b = a.pow(0);
        assert_eq!(b, 1u64);
    }

    #[test]
    fn test_pow_one() {
        let a = Big::from(42u64);
        let b = a.pow(1);
        assert_eq!(b, 42u64);
    }

    #[test]
    fn test_pow_large() {
        let a = Big::from(2u64);
        let b = a.pow(64);
        assert_eq!(format!("{}", b), "18446744073709551616");
    }

    #[test]
    fn test_mod_exp() {
        let base = Big::from(2u64);
        let exp = Big::from(10u64);
        let modulus = Big::from(100u64);
        let result = base.mod_exp(&exp, &modulus);
        assert_eq!(result, 24u64);
    }

    #[test]
    fn test_mod_exp_fermat() {
        let a = Big::from(2u64);
        let p = Big::from(17u64);
        let exp = Big::from(16u64);
        let result = a.mod_exp(&exp, &p);
        assert_eq!(result, 1u64);
    }

    #[test]
    fn test_square() {
        let mut a = Big::from(12u64);
        a.square();
        assert_eq!(a, 144u64);
    }

    #[test]
    fn test_shl_assign() {
        let mut a = Big::from(1u64);
        a.shl_assign(8);
        assert_eq!(a, 256u64);
    }

    #[test]
    fn test_shr_assign() {
        let mut a = Big::from(256u64);
        a.shr_assign(4);
        assert_eq!(a, 16u64);
    }

    #[test]
    fn test_shl_shr_roundtrip() {
        let mut a = Big::from(12345u64);
        a.shl_assign(17);
        a.shr_assign(17);
        assert_eq!(a, 12345u64);
    }

    #[test]
    fn test_large_factorial() {
        let mut result = Big::from(1u64);
        for i in 2..=20u64 {
            result *= i;
        }
        assert_eq!(format!("{}", result), "2432902008176640000");
    }

    #[test]
    fn test_very_large_number() {
        let s = "1234567890".repeat(10);
        let a = Big::from_dec(&s).unwrap();
        assert_eq!(format!("{}", a), s);
    }

    #[test]
    fn test_fibonacci_big() {
        let mut a = Big::from(0u64);
        let mut b = Big::from(1u64);
        for _ in 0..100 {
            let c = &a + &b;
            a = b;
            b = c;
        }
        assert_eq!(format!("{}", a), "354224848179261915075");
    }
}
