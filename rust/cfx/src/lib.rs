//! Rust bindings for cfx - arbitrary precision integers and cryptographic operations.
//!
//! # Example
//! ```
//! use cfx::Big;
//!
//! let a = Big::from(123u64);
//! let b = Big::from(456u64);
//! let c = &a * &b;
//! println!("{} * {} = {}", a, b, c);
//! ```

use std::cmp::Ordering;
use std::ffi::CStr;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign};

/// Raw FFI bindings to the C library.
///
/// These are hand-written for the core functionality.
/// Enable the `bindgen` feature for auto-generated complete bindings.
#[allow(non_camel_case_types)]
pub mod ffi {
    use std::os::raw::{c_char, c_int, c_uint};

    /// The limb type - 64-bit on most platforms.
    #[cfg(target_pointer_width = "64")]
    pub type cfx_limb_t = u64;
    #[cfg(target_pointer_width = "32")]
    pub type cfx_limb_t = u32;

    /// Arbitrary precision positive integer.
    #[repr(C)]
    #[derive(Debug)]
    pub struct cfx_big_t {
        pub limb: *mut cfx_limb_t,
        pub n: usize,
        pub cap: usize,
    }

    extern "C" {
        // Lifecycle
        pub fn cfx_big_init(b: *mut cfx_big_t);
        pub fn cfx_big_free(b: *mut cfx_big_t);
        pub fn cfx_big_copy(dst: *mut cfx_big_t, src: *const cfx_big_t) -> c_int;

        // Conversion from
        pub fn cfx_big_from_u64(b: *mut cfx_big_t, v: u64) -> c_int;
        pub fn cfx_big_from_limb(b: *mut cfx_big_t, v: cfx_limb_t) -> c_int;
        pub fn cfx_big_from_dec(b: *mut cfx_big_t, s: *const c_char) -> c_int;
        pub fn cfx_big_from_hex(b: *mut cfx_big_t, s: *const c_char) -> c_int;
        pub fn cfx_big_from_bytes_be(out: *mut cfx_big_t, be: *const u8, len: usize) -> c_int;

        // Conversion to
        pub fn cfx_big_to_str(b: *const cfx_big_t, sz_out: *mut usize) -> *mut c_char;
        pub fn cfx_big_to_hex(b: *const cfx_big_t, sz_out: *mut usize) -> *mut c_char;

        // Predicates
        pub fn cfx_big_is_zero(b: *const cfx_big_t) -> c_int;
        pub fn cfx_big_is_one(b: *const cfx_big_t) -> c_int;
        pub fn cfx_big_is_even(b: *const cfx_big_t) -> c_int;
        pub fn cfx_big_is_prime(b: *const cfx_big_t) -> c_int;

        // Comparison
        pub fn cfx_big_eq(a: *const cfx_big_t, b: *const cfx_big_t) -> c_int;
        pub fn cfx_big_cmp(a: *const cfx_big_t, b: *const cfx_big_t) -> c_int;
        pub fn cfx_big_cmp_sm(a: *const cfx_big_t, n: cfx_limb_t) -> c_int;

        // Arithmetic (in-place: first arg is modified)
        pub fn cfx_big_add(b: *mut cfx_big_t, a: *const cfx_big_t);
        pub fn cfx_big_add_sm(b: *mut cfx_big_t, n: cfx_limb_t);
        pub fn cfx_big_sub(a: *mut cfx_big_t, b: *const cfx_big_t);
        pub fn cfx_big_sub_sm(b: *mut cfx_big_t, n: cfx_limb_t);
        pub fn cfx_big_mul(out: *mut cfx_big_t, a: *const cfx_big_t, b: *const cfx_big_t);
        pub fn cfx_big_mul_eq(b: *mut cfx_big_t, m: *const cfx_big_t);
        pub fn cfx_big_mul_sm(b: *mut cfx_big_t, m: cfx_limb_t);
        pub fn cfx_big_sq(b: *mut cfx_big_t);
        pub fn cfx_big_div_eq(b: *mut cfx_big_t, d: *const cfx_big_t, r: *mut cfx_big_t) -> c_int;
        pub fn cfx_big_div_sm(b: *mut cfx_big_t, d: cfx_limb_t) -> cfx_limb_t;
        pub fn cfx_big_mod(out: *mut cfx_big_t, n: *const cfx_big_t, m: *const cfx_big_t) -> c_int;
        pub fn cfx_big_mod_sm(b: *const cfx_big_t, m: cfx_limb_t) -> cfx_limb_t;

        // Exponentiation
        pub fn cfx_big_exp(out: *mut cfx_big_t, n: *const cfx_big_t, p: *const cfx_big_t);
        pub fn cfx_big_exp_u64(out: *mut cfx_big_t, n: *const cfx_big_t, p: u64);
        pub fn cfx_big_mod_exp(
            out: *mut cfx_big_t,
            n: *const cfx_big_t,
            p: *const cfx_big_t,
            m: *const cfx_big_t,
        );

        // Bit operations
        pub fn cfx_big_bitlen(b: *const cfx_big_t) -> usize;
        pub fn cfx_big_shl_bits_eq(b: *mut cfx_big_t, s: c_uint);
        pub fn cfx_big_shr_bits_eq(b: *mut cfx_big_t, s: c_uint);
    }
}

/// Arbitrary precision positive integer.
///
/// This is a safe wrapper around `cfx_big_t`.
pub struct Big {
    inner: ffi::cfx_big_t,
}

// Safety: cfx_big_t contains a pointer to heap memory, but doesn't use
// thread-local state. Operations are not thread-safe for the same instance,
// but different instances can be used from different threads.
unsafe impl Send for Big {}

impl Big {
    /// Creates a new Big initialized to zero.
    #[inline]
    pub fn new() -> Self {
        let mut inner = std::mem::MaybeUninit::<ffi::cfx_big_t>::uninit();
        unsafe {
            ffi::cfx_big_init(inner.as_mut_ptr());
            Big {
                inner: inner.assume_init(),
            }
        }
    }

    /// Returns true if this number is zero.
    #[inline]
    pub fn is_zero(&self) -> bool {
        unsafe { ffi::cfx_big_is_zero(&self.inner) != 0 }
    }

    /// Returns true if this number is one.
    #[inline]
    pub fn is_one(&self) -> bool {
        unsafe { ffi::cfx_big_is_one(&self.inner) != 0 }
    }

    /// Returns true if this number is even.
    #[inline]
    pub fn is_even(&self) -> bool {
        unsafe { ffi::cfx_big_is_even(&self.inner) != 0 }
    }

    /// Returns true if this number is odd.
    #[inline]
    pub fn is_odd(&self) -> bool {
        !self.is_even()
    }

    /// Returns true if this number is probably prime (Miller-Rabin).
    #[inline]
    pub fn is_prime(&self) -> bool {
        unsafe { ffi::cfx_big_is_prime(&self.inner) != 0 }
    }

    /// Returns the number of bits needed to represent this number.
    #[inline]
    pub fn bit_len(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            unsafe { ffi::cfx_big_bitlen(&self.inner) }
        }
    }

    /// Computes `self^exp mod modulus`.
    pub fn mod_exp(&self, exp: &Big, modulus: &Big) -> Big {
        let mut result = Big::new();
        unsafe {
            ffi::cfx_big_mod_exp(&mut result.inner, &self.inner, &exp.inner, &modulus.inner);
        }
        result
    }

    /// Computes `self^exp`.
    pub fn pow(&self, exp: u64) -> Big {
        let mut result = Big::new();
        unsafe {
            ffi::cfx_big_exp_u64(&mut result.inner, &self.inner, exp);
        }
        result
    }

    /// Computes `self * self` (squaring).
    pub fn square(&mut self) {
        unsafe {
            ffi::cfx_big_sq(&mut self.inner);
        }
    }

    /// Left shift by `bits` bits (multiply by 2^bits).
    pub fn shl_assign(&mut self, bits: u32) {
        unsafe {
            ffi::cfx_big_shl_bits_eq(&mut self.inner, bits);
        }
    }

    /// Right shift by `bits` bits (divide by 2^bits).
    pub fn shr_assign(&mut self, bits: u32) {
        unsafe {
            ffi::cfx_big_shr_bits_eq(&mut self.inner, bits);
        }
    }

    /// Parses a decimal string into a Big.
    pub fn from_dec(s: &str) -> Option<Big> {
        let c_str = std::ffi::CString::new(s).ok()?;
        let mut result = Big::new();
        let rc = unsafe { ffi::cfx_big_from_dec(&mut result.inner, c_str.as_ptr()) };
        if rc == 0 {
            Some(result)
        } else {
            None
        }
    }

    /// Parses a hexadecimal string into a Big.
    pub fn from_hex(s: &str) -> Option<Big> {
        let c_str = std::ffi::CString::new(s).ok()?;
        let mut result = Big::new();
        let rc = unsafe { ffi::cfx_big_from_hex(&mut result.inner, c_str.as_ptr()) };
        if rc == 0 {
            Some(result)
        } else {
            None
        }
    }

    /// Creates a Big from big-endian bytes.
    pub fn from_bytes_be(bytes: &[u8]) -> Option<Big> {
        let mut result = Big::new();
        let rc = unsafe {
            ffi::cfx_big_from_bytes_be(&mut result.inner, bytes.as_ptr(), bytes.len())
        };
        if rc == 0 {
            Some(result)
        } else {
            None
        }
    }

    /// Returns `self % m` where m is a u64.
    #[inline]
    pub fn mod_sm(&self, m: u64) -> u64 {
        unsafe { ffi::cfx_big_mod_sm(&self.inner, m as ffi::cfx_limb_t) as u64 }
    }

    /// In-place addition of a u64: `self += n`.
    #[inline]
    pub fn add_sm(&mut self, n: u64) {
        unsafe { ffi::cfx_big_add_sm(&mut self.inner, n as ffi::cfx_limb_t) };
    }

    /// In-place subtraction of a u64: `self -= n`.
    #[inline]
    pub fn sub_sm(&mut self, n: u64) {
        unsafe { ffi::cfx_big_sub_sm(&mut self.inner, n as ffi::cfx_limb_t) };
    }

    /// Returns the hexadecimal string representation.
    pub fn to_hex(&self) -> String {
        let mut len: usize = 0;
        let ptr = unsafe { ffi::cfx_big_to_hex(&self.inner, &mut len) };
        if ptr.is_null() {
            return String::new();
        }
        let s = unsafe { CStr::from_ptr(ptr).to_string_lossy().into_owned() };
        unsafe { libc_free(ptr as *mut _) };
        s
    }

    /// Provides access to the underlying C type (unsafe).
    #[inline]
    pub unsafe fn as_ptr(&self) -> *const ffi::cfx_big_t {
        &self.inner
    }

    /// Provides mutable access to the underlying C type (unsafe).
    #[inline]
    pub unsafe fn as_mut_ptr(&mut self) -> *mut ffi::cfx_big_t {
        &mut self.inner
    }
}

// Free memory allocated by C
#[cfg(unix)]
unsafe fn libc_free(ptr: *mut std::ffi::c_void) {
    extern "C" {
        fn free(ptr: *mut std::ffi::c_void);
    }
    free(ptr);
}

#[cfg(windows)]
unsafe fn libc_free(ptr: *mut std::ffi::c_void) {
    extern "C" {
        fn free(ptr: *mut std::ffi::c_void);
    }
    free(ptr);
}

impl Default for Big {
    fn default() -> Self {
        Big::new()
    }
}

impl Drop for Big {
    fn drop(&mut self) {
        unsafe {
            ffi::cfx_big_free(&mut self.inner);
        }
    }
}

impl Clone for Big {
    fn clone(&self) -> Self {
        let mut new = Big::new();
        unsafe {
            ffi::cfx_big_copy(&mut new.inner, &self.inner);
        }
        new
    }
}

impl fmt::Debug for Big {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Big({})", self)
    }
}

impl fmt::Display for Big {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut len: usize = 0;
        let ptr = unsafe { ffi::cfx_big_to_str(&self.inner, &mut len) };
        if ptr.is_null() {
            return write!(f, "0");
        }
        let s = unsafe { CStr::from_ptr(ptr).to_string_lossy() };
        let result = write!(f, "{}", s);
        unsafe { libc_free(ptr as *mut _) };
        result
    }
}

// ============ From implementations ============

impl From<u64> for Big {
    fn from(v: u64) -> Self {
        let mut b = Big::new();
        unsafe {
            ffi::cfx_big_from_u64(&mut b.inner, v);
        }
        b
    }
}

impl From<u32> for Big {
    fn from(v: u32) -> Self {
        Big::from(v as u64)
    }
}

impl From<usize> for Big {
    fn from(v: usize) -> Self {
        Big::from(v as u64)
    }
}

// ============ Comparison ============

impl PartialEq for Big {
    fn eq(&self, other: &Self) -> bool {
        unsafe { ffi::cfx_big_eq(&self.inner, &other.inner) != 0 }
    }
}

impl Eq for Big {}

impl PartialOrd for Big {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Big {
    fn cmp(&self, other: &Self) -> Ordering {
        let r = unsafe { ffi::cfx_big_cmp(&self.inner, &other.inner) };
        match r {
            -1 => Ordering::Less,
            0 => Ordering::Equal,
            1 => Ordering::Greater,
            _ => Ordering::Equal, // shouldn't happen
        }
    }
}

impl PartialEq<u64> for Big {
    fn eq(&self, other: &u64) -> bool {
        unsafe { ffi::cfx_big_cmp_sm(&self.inner, *other as ffi::cfx_limb_t) == 0 }
    }
}

// ============ Arithmetic (consuming) ============

impl Add for Big {
    type Output = Big;
    fn add(mut self, rhs: Big) -> Big {
        unsafe { ffi::cfx_big_add(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a> Add<&'a Big> for Big {
    type Output = Big;
    fn add(mut self, rhs: &'a Big) -> Big {
        unsafe { ffi::cfx_big_add(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a, 'b> Add<&'b Big> for &'a Big {
    type Output = Big;
    fn add(self, rhs: &'b Big) -> Big {
        let mut result = self.clone();
        unsafe { ffi::cfx_big_add(&mut result.inner, &rhs.inner) };
        result
    }
}

impl Add<u64> for Big {
    type Output = Big;
    fn add(mut self, rhs: u64) -> Big {
        unsafe { ffi::cfx_big_add_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
        self
    }
}

impl Sub for Big {
    type Output = Big;
    fn sub(mut self, rhs: Big) -> Big {
        unsafe { ffi::cfx_big_sub(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a> Sub<&'a Big> for Big {
    type Output = Big;
    fn sub(mut self, rhs: &'a Big) -> Big {
        unsafe { ffi::cfx_big_sub(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a, 'b> Sub<&'b Big> for &'a Big {
    type Output = Big;
    fn sub(self, rhs: &'b Big) -> Big {
        let mut result = self.clone();
        unsafe { ffi::cfx_big_sub(&mut result.inner, &rhs.inner) };
        result
    }
}

impl Sub<u64> for Big {
    type Output = Big;
    fn sub(mut self, rhs: u64) -> Big {
        unsafe { ffi::cfx_big_sub_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
        self
    }
}

impl Mul for Big {
    type Output = Big;
    fn mul(mut self, rhs: Big) -> Big {
        unsafe { ffi::cfx_big_mul_eq(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a> Mul<&'a Big> for Big {
    type Output = Big;
    fn mul(mut self, rhs: &'a Big) -> Big {
        unsafe { ffi::cfx_big_mul_eq(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a, 'b> Mul<&'b Big> for &'a Big {
    type Output = Big;
    fn mul(self, rhs: &'b Big) -> Big {
        let mut result = self.clone();
        unsafe { ffi::cfx_big_mul_eq(&mut result.inner, &rhs.inner) };
        result
    }
}

impl Mul<u64> for Big {
    type Output = Big;
    fn mul(mut self, rhs: u64) -> Big {
        unsafe { ffi::cfx_big_mul_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
        self
    }
}

impl Div for Big {
    type Output = Big;
    fn div(mut self, rhs: Big) -> Big {
        unsafe { ffi::cfx_big_div_eq(&mut self.inner, &rhs.inner, std::ptr::null_mut()) };
        self
    }
}

impl<'a> Div<&'a Big> for Big {
    type Output = Big;
    fn div(mut self, rhs: &'a Big) -> Big {
        unsafe { ffi::cfx_big_div_eq(&mut self.inner, &rhs.inner, std::ptr::null_mut()) };
        self
    }
}

impl<'a, 'b> Div<&'b Big> for &'a Big {
    type Output = Big;
    fn div(self, rhs: &'b Big) -> Big {
        let mut result = self.clone();
        unsafe { ffi::cfx_big_div_eq(&mut result.inner, &rhs.inner, std::ptr::null_mut()) };
        result
    }
}

impl Div<u64> for Big {
    type Output = Big;
    fn div(mut self, rhs: u64) -> Big {
        unsafe { ffi::cfx_big_div_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
        self
    }
}

impl Rem for Big {
    type Output = Big;
    fn rem(self, rhs: Big) -> Big {
        let mut result = Big::new();
        unsafe { ffi::cfx_big_mod(&mut result.inner, &self.inner, &rhs.inner) };
        result
    }
}

impl<'a> Rem<&'a Big> for Big {
    type Output = Big;
    fn rem(self, rhs: &'a Big) -> Big {
        let mut result = Big::new();
        unsafe { ffi::cfx_big_mod(&mut result.inner, &self.inner, &rhs.inner) };
        result
    }
}

impl<'a, 'b> Rem<&'b Big> for &'a Big {
    type Output = Big;
    fn rem(self, rhs: &'b Big) -> Big {
        let mut result = Big::new();
        unsafe { ffi::cfx_big_mod(&mut result.inner, &self.inner, &rhs.inner) };
        result
    }
}

impl Rem<u64> for Big {
    type Output = u64;
    fn rem(self, rhs: u64) -> u64 {
        unsafe { ffi::cfx_big_mod_sm(&self.inner, rhs as ffi::cfx_limb_t) as u64 }
    }
}

impl<'a> Rem<u64> for &'a Big {
    type Output = u64;
    fn rem(self, rhs: u64) -> u64 {
        unsafe { ffi::cfx_big_mod_sm(&self.inner, rhs as ffi::cfx_limb_t) as u64 }
    }
}

// ============ Assign variants ============

impl AddAssign<&Big> for Big {
    fn add_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_add(&mut self.inner, &rhs.inner) };
    }
}

impl AddAssign<u64> for Big {
    fn add_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_add_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
    }
}

impl SubAssign<&Big> for Big {
    fn sub_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_sub(&mut self.inner, &rhs.inner) };
    }
}

impl SubAssign<u64> for Big {
    fn sub_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_sub_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
    }
}

impl MulAssign<&Big> for Big {
    fn mul_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_mul_eq(&mut self.inner, &rhs.inner) };
    }
}

impl MulAssign<u64> for Big {
    fn mul_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_mul_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
    }
}

impl DivAssign<&Big> for Big {
    fn div_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_div_eq(&mut self.inner, &rhs.inner, std::ptr::null_mut()) };
    }
}

impl DivAssign<u64> for Big {
    fn div_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_div_sm(&mut self.inner, rhs as ffi::cfx_limb_t) };
    }
}

impl RemAssign<&Big> for Big {
    fn rem_assign(&mut self, rhs: &Big) {
        let mut result = Big::new();
        unsafe { ffi::cfx_big_mod(&mut result.inner, &self.inner, &rhs.inner) };
        std::mem::swap(self, &mut result);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ============ Basic construction and display ============

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

    // ============ String parsing ============

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
        assert!(Big::from_dec("123abc").is_none());
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

    // ============ Predicates ============

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
        // Larger primes
        assert!(Big::from(104729u64).is_prime()); // 10000th prime
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

    // ============ Comparison ============

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

    // ============ Addition ============

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

    // ============ Subtraction ============

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

    // ============ Multiplication ============

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
        // Test multiplication with large numbers
        let a = Big::from_dec("123456789012345678901234567890").unwrap();
        let b = Big::from(2u64);
        let c = &a * &b;
        assert_eq!(format!("{}", c), "246913578024691357802469135780");
    }

    // ============ Division ============

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

    // ============ Remainder ============

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

    // ============ Power and modular exponentiation ============

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
        assert_eq!(format!("{}", b), "18446744073709551616"); // 2^64
    }

    #[test]
    fn test_mod_exp() {
        // 2^10 mod 100 = 1024 mod 100 = 24
        let base = Big::from(2u64);
        let exp = Big::from(10u64);
        let modulus = Big::from(100u64);
        let result = base.mod_exp(&exp, &modulus);
        assert_eq!(result, 24u64);
    }

    #[test]
    fn test_mod_exp_fermat() {
        // Fermat's little theorem: a^(p-1) ≡ 1 (mod p) for prime p
        let a = Big::from(2u64);
        let p = Big::from(17u64); // prime
        let exp = Big::from(16u64); // p - 1
        let result = a.mod_exp(&exp, &p);
        assert_eq!(result, 1u64);
    }

    // ============ Square and shift ============

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

    // ============ Edge cases and large numbers ============

    #[test]
    fn test_large_factorial() {
        // Compute 20!
        let mut result = Big::from(1u64);
        for i in 2..=20u64 {
            result *= i;
        }
        assert_eq!(format!("{}", result), "2432902008176640000");
    }

    #[test]
    fn test_very_large_number() {
        // Create a 100-digit number
        let s = "1234567890".repeat(10);
        let a = Big::from_dec(&s).unwrap();
        assert_eq!(format!("{}", a), s);
    }

    #[test]
    fn test_fibonacci_big() {
        // Compute F(100) - a 21-digit number
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
