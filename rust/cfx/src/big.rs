//! Arbitrary precision positive integers.
//!
//! This module provides the `Big` type, a safe wrapper around the C library's
//! `cfx_big_t` for arbitrary precision arithmetic.
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

use crate::ffi;
use std::cmp::Ordering;
use std::ffi::CStr;
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign};

/// Arbitrary precision positive integer.
///
/// This is a safe wrapper around `cfx_big_t`. Values are always non-negative.
pub struct Big {
    pub(crate) inner: ffi::cfx_big_t,
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
            ffi::cfx_big_sq_eq(&mut self.inner);
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
        unsafe { ffi::cfx_big_add_sm_eq(&mut self.inner, n as ffi::cfx_limb_t) };
    }

    /// In-place subtraction of a u64: `self -= n`.
    #[inline]
    pub fn sub_sm(&mut self, n: u64) {
        unsafe { ffi::cfx_big_sub_sm_eq(&mut self.inner, n as ffi::cfx_limb_t) };
    }

    /// Returns the hexadecimal string representation.
    pub fn to_hex(&self) -> String {
        let mut len: usize = 0;
        let ptr = unsafe { ffi::cfx_big_hex_alloc(&self.inner, &mut len) };
        if ptr.is_null() {
            return String::new();
        }
        let s = unsafe { CStr::from_ptr(ptr).to_string_lossy().into_owned() };
        unsafe { ffi::libc_free(ptr as *mut _) };
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
        let ptr = unsafe { ffi::cfx_big_dec_alloc(&self.inner, &mut len) };
        if ptr.is_null() {
            return write!(f, "0");
        }
        let s = unsafe { CStr::from_ptr(ptr).to_string_lossy() };
        let result = write!(f, "{}", s);
        unsafe { ffi::libc_free(ptr as *mut _) };
        result
    }
}

// From implementations

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

// Comparison

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
            _ => Ordering::Equal,
        }
    }
}

impl PartialEq<u64> for Big {
    fn eq(&self, other: &u64) -> bool {
        unsafe { ffi::cfx_big_cmp_sm(&self.inner, *other as ffi::cfx_limb_t) == 0 }
    }
}

// Arithmetic (consuming)

impl Add for Big {
    type Output = Big;
    fn add(mut self, rhs: Big) -> Big {
        unsafe { ffi::cfx_big_add_eq(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a> Add<&'a Big> for Big {
    type Output = Big;
    fn add(mut self, rhs: &'a Big) -> Big {
        unsafe { ffi::cfx_big_add_eq(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a, 'b> Add<&'b Big> for &'a Big {
    type Output = Big;
    fn add(self, rhs: &'b Big) -> Big {
        let mut result = self.clone();
        unsafe { ffi::cfx_big_add_eq(&mut result.inner, &rhs.inner) };
        result
    }
}

impl Add<u64> for Big {
    type Output = Big;
    fn add(mut self, rhs: u64) -> Big {
        unsafe { ffi::cfx_big_add_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
        self
    }
}

impl Sub for Big {
    type Output = Big;
    fn sub(mut self, rhs: Big) -> Big {
        unsafe { ffi::cfx_big_sub_eq(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a> Sub<&'a Big> for Big {
    type Output = Big;
    fn sub(mut self, rhs: &'a Big) -> Big {
        unsafe { ffi::cfx_big_sub_eq(&mut self.inner, &rhs.inner) };
        self
    }
}

impl<'a, 'b> Sub<&'b Big> for &'a Big {
    type Output = Big;
    fn sub(self, rhs: &'b Big) -> Big {
        let mut result = self.clone();
        unsafe { ffi::cfx_big_sub_eq(&mut result.inner, &rhs.inner) };
        result
    }
}

impl Sub<u64> for Big {
    type Output = Big;
    fn sub(mut self, rhs: u64) -> Big {
        unsafe { ffi::cfx_big_sub_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
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
        unsafe { ffi::cfx_big_mul_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
        self
    }
}

impl Div for Big {
    type Output = Big;
    fn div(mut self, rhs: Big) -> Big {
        unsafe { ffi::cfx_big_divrem_eq(&mut self.inner, &rhs.inner, std::ptr::null_mut()) };
        self
    }
}

impl<'a> Div<&'a Big> for Big {
    type Output = Big;
    fn div(mut self, rhs: &'a Big) -> Big {
        unsafe { ffi::cfx_big_divrem_eq(&mut self.inner, &rhs.inner, std::ptr::null_mut()) };
        self
    }
}

impl<'a, 'b> Div<&'b Big> for &'a Big {
    type Output = Big;
    fn div(self, rhs: &'b Big) -> Big {
        let mut result = self.clone();
        unsafe { ffi::cfx_big_divrem_eq(&mut result.inner, &rhs.inner, std::ptr::null_mut()) };
        result
    }
}

impl Div<u64> for Big {
    type Output = Big;
    fn div(mut self, rhs: u64) -> Big {
        unsafe { ffi::cfx_big_div_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
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

// Assign variants

impl AddAssign<&Big> for Big {
    fn add_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_add_eq(&mut self.inner, &rhs.inner) };
    }
}

impl AddAssign<u64> for Big {
    fn add_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_add_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
    }
}

impl SubAssign<&Big> for Big {
    fn sub_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_sub_eq(&mut self.inner, &rhs.inner) };
    }
}

impl SubAssign<u64> for Big {
    fn sub_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_sub_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
    }
}

impl MulAssign<&Big> for Big {
    fn mul_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_mul_eq(&mut self.inner, &rhs.inner) };
    }
}

impl MulAssign<u64> for Big {
    fn mul_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_mul_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
    }
}

impl DivAssign<&Big> for Big {
    fn div_assign(&mut self, rhs: &Big) {
        unsafe { ffi::cfx_big_divrem_eq(&mut self.inner, &rhs.inner, std::ptr::null_mut()) };
    }
}

impl DivAssign<u64> for Big {
    fn div_assign(&mut self, rhs: u64) {
        unsafe { ffi::cfx_big_div_sm_eq(&mut self.inner, rhs as ffi::cfx_limb_t) };
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
    fn test_arithmetic() {
        let a = Big::from(100u64);
        let b = Big::from(23u64);
        assert_eq!(&a + &b, 123u64);
        assert_eq!(&a - &b, 77u64);
        assert_eq!(&a * &b, 2300u64);
    }

    #[test]
    fn test_pow() {
        let a = Big::from(2u64);
        assert_eq!(a.pow(10), 1024u64);
    }

    #[test]
    fn test_mod_exp() {
        let base = Big::from(2u64);
        let exp = Big::from(10u64);
        let modulus = Big::from(100u64);
        assert_eq!(base.mod_exp(&exp, &modulus), 24u64);
    }
}
