//! Safe, single-threaded ownership wrapper around unmodified libgap.

use crate::Error;
use std::cell::RefCell;
use std::env;
use std::ffi::{c_char, CStr, CString};
use std::marker::PhantomData;
use std::path::Path;
use std::rc::Rc;
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_SLOT: AtomicUsize = AtomicUsize::new(1);

mod ffi {
    use super::c_char;

    unsafe extern "C" {
        pub fn libgap_cyc_init(gap_root: *const c_char) -> i32;
        pub fn libgap_cyc_error() -> *const c_char;
        pub fn libgap_cyc_from_terms(
            slot: usize,
            order: u32,
            len: usize,
            exponents: *const u32,
            numerators: *const i64,
            denominators: *const i64,
        ) -> i32;
        pub fn libgap_cyc_add(output: usize, lhs: usize, rhs: usize) -> i32;
        pub fn libgap_cyc_mul(output: usize, lhs: usize, rhs: usize) -> i32;
        pub fn libgap_cyc_quo(output: usize, lhs: usize, rhs: usize) -> i32;
        pub fn libgap_cyc_release(slot: usize) -> i32;
        pub fn libgap_cyc_equal(lhs: usize, rhs: usize, equal: *mut i32) -> i32;
        pub fn libgap_cyc_order(slot: usize, order: *mut u32) -> i32;
        pub fn libgap_cyc_coefficient(
            slot: usize,
            exponent: u32,
            numerator: *mut i64,
            denominator: *mut i64,
        ) -> i32;
    }
}

fn error() -> Error {
    let pointer = unsafe { ffi::libgap_cyc_error() };
    if pointer.is_null() {
        return Error("unknown libgap error".into());
    }
    Error(
        unsafe { CStr::from_ptr(pointer) }
            .to_string_lossy()
            .into_owned(),
    )
}

struct Slots {
    free: Vec<usize>,
}

impl Slots {
    fn allocate(&mut self) -> usize {
        self.free
            .pop()
            .unwrap_or_else(|| NEXT_SLOT.fetch_add(1, Ordering::Relaxed))
    }
}

pub struct Context {
    slots: Rc<RefCell<Slots>>,
    _not_send_or_sync: PhantomData<Rc<()>>,
}

pub struct Cyclotomic {
    slot: usize,
    slots: Rc<RefCell<Slots>>,
    _not_send_or_sync: PhantomData<Rc<()>>,
}

impl Context {
    pub fn new() -> Result<Self, Error> {
        let root =
            env::var_os("LIBGAP_ROOT").ok_or_else(|| Error("LIBGAP_ROOT is not set".into()))?;
        Self::new_at(root)
    }

    pub fn new_at(root: impl AsRef<Path>) -> Result<Self, Error> {
        let root = CString::new(root.as_ref().as_os_str().as_encoded_bytes())
            .map_err(|_| Error("GAP root contains an interior NUL byte".into()))?;
        if unsafe { ffi::libgap_cyc_init(root.as_ptr()) } == 0 {
            return Err(error());
        }
        Ok(Self {
            slots: Rc::new(RefCell::new(Slots { free: Vec::new() })),
            _not_send_or_sync: PhantomData,
        })
    }

    fn output(&self, operation: impl FnOnce(usize) -> i32) -> Result<Cyclotomic, Error> {
        let slot = self.slots.borrow_mut().allocate();
        if operation(slot) == 0 {
            self.slots.borrow_mut().free.push(slot);
            return Err(error());
        }
        Ok(Cyclotomic {
            slot,
            slots: Rc::clone(&self.slots),
            _not_send_or_sync: PhantomData,
        })
    }

    pub fn from_terms(&self, order: u32, terms: &[(u32, (i64, i64))]) -> Result<Cyclotomic, Error> {
        if order == 0 {
            return Err(Error("cyclotomic order must be positive".into()));
        }
        let exponents: Vec<_> = terms.iter().map(|term| term.0).collect();
        let numerators: Vec<_> = terms.iter().map(|term| term.1 .0).collect();
        let denominators: Vec<_> = terms.iter().map(|term| term.1 .1).collect();
        self.output(|slot| unsafe {
            ffi::libgap_cyc_from_terms(
                slot,
                order,
                terms.len(),
                exponents.as_ptr(),
                numerators.as_ptr(),
                denominators.as_ptr(),
            )
        })
    }

    pub fn from_integer_terms(
        &self,
        order: u32,
        terms: &[(u32, i64)],
    ) -> Result<Cyclotomic, Error> {
        let rational_terms: Vec<_> = terms
            .iter()
            .map(|&(exponent, coefficient)| (exponent, (coefficient, 1)))
            .collect();
        self.from_terms(order, &rational_terms)
    }

    pub fn add(&self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        self.output(|slot| unsafe { ffi::libgap_cyc_add(slot, lhs.slot, rhs.slot) })
    }

    pub fn mul(&self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        self.output(|slot| unsafe { ffi::libgap_cyc_mul(slot, lhs.slot, rhs.slot) })
    }

    pub fn quo(&self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        self.output(|slot| unsafe { ffi::libgap_cyc_quo(slot, lhs.slot, rhs.slot) })
    }
}

impl Cyclotomic {
    pub fn order(&self) -> Result<u32, Error> {
        let mut order = 0;
        if unsafe { ffi::libgap_cyc_order(self.slot, &mut order) } == 0 {
            return Err(error());
        }
        Ok(order)
    }

    pub fn coefficients(&self) -> Result<Vec<(i64, i64)>, Error> {
        let order = self.order()?;
        (0..order)
            .map(|exponent| {
                let mut numerator = 0;
                let mut denominator = 0;
                if unsafe {
                    ffi::libgap_cyc_coefficient(
                        self.slot,
                        exponent,
                        &mut numerator,
                        &mut denominator,
                    )
                } == 0
                {
                    return Err(error());
                }
                Ok((numerator, denominator))
            })
            .collect()
    }
}

impl Drop for Cyclotomic {
    fn drop(&mut self) {
        unsafe {
            ffi::libgap_cyc_release(self.slot);
        }
        self.slots.borrow_mut().free.push(self.slot);
    }
}

impl PartialEq for Cyclotomic {
    fn eq(&self, other: &Self) -> bool {
        let mut equal = 0;
        unsafe { ffi::libgap_cyc_equal(self.slot, other.slot, &mut equal) != 0 && equal != 0 }
    }
}

impl Eq for Cyclotomic {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unmodified_gap_handles_rational_cyclotomics() {
        let context = Context::new().unwrap();
        let value = context.from_terms(5, &[(1, (1, 2)), (2, (2, 3))]).unwrap();
        let square = context.mul(&value, &value).unwrap();
        assert_eq!(square.order().unwrap(), 5);
        assert_eq!(
            square.coefficients().unwrap(),
            vec![(0, 1), (0, 1), (1, 4), (2, 3), (4, 9)]
        );
        let recovered = context.quo(&square, &value).unwrap();
        assert!(recovered == value);
    }
}
