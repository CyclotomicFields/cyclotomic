//! Safe Rust ownership layer for the standalone GAP-derived C reference.
//!
//! This crate intentionally supports only checked `i64` coefficients in its
//! first milestone. It is a benchmark and differential-testing reference, not
//! a production API.

use std::ffi::{c_char, CStr};
use std::fmt;
use std::marker::PhantomData;
use std::ptr::NonNull;
use std::rc::Rc;

#[cfg(feature = "libgap")]
pub mod libgap;

pub mod literal;
pub mod tagged;
#[cfg(feature = "libgap")]
pub mod tagged_character;

mod ffi {
    use super::c_char;

    #[repr(C)]
    pub struct GapCycCtx {
        _private: [u8; 0],
    }

    #[repr(C)]
    pub struct GapCyc {
        _private: [u8; 0],
    }

    unsafe extern "C" {
        pub fn gap_cyc_ctx_new() -> *mut GapCycCtx;
        pub fn gap_cyc_ctx_free(ctx: *mut GapCycCtx);
        pub fn gap_cyc_ctx_error(ctx: *const GapCycCtx) -> *const c_char;

        pub fn gap_cyc_from_terms(
            ctx: *mut GapCycCtx,
            order: u32,
            len: usize,
            exponents: *const u32,
            coefficients: *const i64,
        ) -> *mut GapCyc;
        pub fn gap_cyc_root(ctx: *mut GapCycCtx, order: u32, exponent: u32) -> *mut GapCyc;
        pub fn gap_cyc_add(
            ctx: *mut GapCycCtx,
            lhs: *const GapCyc,
            rhs: *const GapCyc,
        ) -> *mut GapCyc;
        pub fn gap_cyc_mul(
            ctx: *mut GapCycCtx,
            lhs: *const GapCyc,
            rhs: *const GapCyc,
        ) -> *mut GapCyc;

        pub fn gap_cyc_free(cyc: *mut GapCyc);
        pub fn gap_cyc_order(cyc: *const GapCyc) -> u32;
        pub fn gap_cyc_len(cyc: *const GapCyc) -> usize;
        pub fn gap_cyc_exponent(cyc: *const GapCyc, index: usize) -> u32;
        pub fn gap_cyc_coefficient(cyc: *const GapCyc, index: usize) -> i64;
        pub fn gap_cyc_equal(lhs: *const GapCyc, rhs: *const GapCyc) -> i32;
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Error(String);

impl fmt::Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

pub struct Context {
    raw: NonNull<ffi::GapCycCtx>,
    _not_send_or_sync: PhantomData<Rc<()>>,
}

pub struct Cyclotomic {
    raw: NonNull<ffi::GapCyc>,
    _not_send_or_sync: PhantomData<Rc<()>>,
}

impl Context {
    pub fn new() -> Result<Self, Error> {
        let raw = NonNull::new(unsafe { ffi::gap_cyc_ctx_new() })
            .ok_or_else(|| Error("failed to allocate GAP cyclotomic context".into()))?;
        Ok(Self {
            raw,
            _not_send_or_sync: PhantomData,
        })
    }

    fn error(&self) -> Error {
        let pointer = unsafe { ffi::gap_cyc_ctx_error(self.raw.as_ptr()) };
        if pointer.is_null() {
            return Error("unknown GAP cyclotomic error".into());
        }
        let message = unsafe { CStr::from_ptr(pointer) }
            .to_string_lossy()
            .into_owned();
        Error(message)
    }

    fn element(&self, raw: *mut ffi::GapCyc) -> Result<Cyclotomic, Error> {
        let raw = NonNull::new(raw).ok_or_else(|| self.error())?;
        Ok(Cyclotomic {
            raw,
            _not_send_or_sync: PhantomData,
        })
    }

    pub fn from_terms(&mut self, order: u32, terms: &[(u32, i64)]) -> Result<Cyclotomic, Error> {
        let exponents: Vec<_> = terms.iter().map(|term| term.0).collect();
        let coefficients: Vec<_> = terms.iter().map(|term| term.1).collect();
        let raw = unsafe {
            ffi::gap_cyc_from_terms(
                self.raw.as_ptr(),
                order,
                terms.len(),
                exponents.as_ptr(),
                coefficients.as_ptr(),
            )
        };
        self.element(raw)
    }

    pub fn root(&mut self, order: u32, exponent: u32) -> Result<Cyclotomic, Error> {
        let raw = unsafe { ffi::gap_cyc_root(self.raw.as_ptr(), order, exponent) };
        self.element(raw)
    }

    pub fn add(&mut self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        let raw =
            unsafe { ffi::gap_cyc_add(self.raw.as_ptr(), lhs.raw.as_ptr(), rhs.raw.as_ptr()) };
        self.element(raw)
    }

    pub fn mul(&mut self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        let raw =
            unsafe { ffi::gap_cyc_mul(self.raw.as_ptr(), lhs.raw.as_ptr(), rhs.raw.as_ptr()) };
        self.element(raw)
    }
}

impl Drop for Context {
    fn drop(&mut self) {
        unsafe { ffi::gap_cyc_ctx_free(self.raw.as_ptr()) }
    }
}

impl Cyclotomic {
    pub fn order(&self) -> u32 {
        unsafe { ffi::gap_cyc_order(self.raw.as_ptr()) }
    }

    pub fn terms(&self) -> Vec<(u32, i64)> {
        let len = unsafe { ffi::gap_cyc_len(self.raw.as_ptr()) };
        (0..len)
            .map(|index| unsafe {
                (
                    ffi::gap_cyc_exponent(self.raw.as_ptr(), index),
                    ffi::gap_cyc_coefficient(self.raw.as_ptr(), index),
                )
            })
            .collect()
    }
}

impl Drop for Cyclotomic {
    fn drop(&mut self) {
        unsafe { ffi::gap_cyc_free(self.raw.as_ptr()) }
    }
}

impl PartialEq for Cyclotomic {
    fn eq(&self, other: &Self) -> bool {
        unsafe { ffi::gap_cyc_equal(self.raw.as_ptr(), other.raw.as_ptr()) != 0 }
    }
}

impl Eq for Cyclotomic {}

impl fmt::Debug for Cyclotomic {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Cyclotomic")
            .field("order", &self.order())
            .field("terms", &self.terms())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn third_roots_sum_to_zero() {
        let mut context = Context::new().unwrap();
        let one = context.root(3, 0).unwrap();
        let e = context.root(3, 1).unwrap();
        let e2 = context.root(3, 2).unwrap();
        let one_plus_e = context.add(&one, &e).unwrap();
        let sum = context.add(&one_plus_e, &e2).unwrap();
        assert_eq!(sum.order(), 1);
        assert!(sum.terms().is_empty());
    }

    #[test]
    fn multiplication_reduces_exponents() {
        let mut context = Context::new().unwrap();
        let e2 = context.root(5, 2).unwrap();
        let e4 = context.root(5, 4).unwrap();
        let product = context.mul(&e2, &e4).unwrap();
        let e = context.root(5, 1).unwrap();
        assert_eq!(product, e);
    }

    #[test]
    fn reduces_to_a_smaller_conductor() {
        let mut context = Context::new().unwrap();
        let e6 = context.root(6, 1).unwrap();
        assert_eq!(e6.order(), 3);
    }

    #[test]
    fn combines_different_fields() {
        let mut context = Context::new().unwrap();
        let e3 = context.root(3, 1).unwrap();
        let e5 = context.root(5, 1).unwrap();
        let product = context.mul(&e3, &e5).unwrap();
        assert_eq!(product.order(), 15);
    }
}
