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
        pub fn libgap_character_table_dimensions(
            table: u32,
            rows: *mut usize,
            columns: *mut usize,
        ) -> i32;
        pub fn libgap_character_table(
            table: u32,
            slots: *const usize,
            slots_len: usize,
            class_sizes: *mut i64,
            group_order: *mut i64,
        ) -> i32;
        pub fn libgap_character_tensor_decomposition(
            table: u32,
            lhs: usize,
            rhs: usize,
            multiplicities: *mut i64,
            multiplicities_len: usize,
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

#[derive(Clone, Copy, Debug)]
#[repr(u32)]
pub enum CharacterTableCase {
    A5 = 0,
    Sl2_5 = 1,
    Psl2_11 = 2,
}

impl CharacterTableCase {
    pub const ALL: [Self; 3] = [Self::A5, Self::Sl2_5, Self::Psl2_11];

    pub fn name(self) -> &'static str {
        match self {
            Self::A5 => "A5",
            Self::Sl2_5 => "SL(2,5)",
            Self::Psl2_11 => "PSL(2,11)",
        }
    }
}

pub struct CharacterTable {
    pub case: CharacterTableCase,
    pub rows: Vec<Vec<Cyclotomic>>,
    pub class_sizes: Vec<i64>,
    pub group_order: i64,
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

    pub fn character_table(&self, case: CharacterTableCase) -> Result<CharacterTable, Error> {
        let mut row_count = 0;
        let mut column_count = 0;
        if unsafe {
            ffi::libgap_character_table_dimensions(case as u32, &mut row_count, &mut column_count)
        } == 0
        {
            return Err(error());
        }

        let slot_count = row_count
            .checked_mul(column_count)
            .ok_or_else(|| Error("character table is too large".into()))?;
        let slots: Vec<_> = (0..slot_count)
            .map(|_| self.slots.borrow_mut().allocate())
            .collect();
        let mut class_sizes = vec![0; column_count];
        let mut group_order = 0;
        if unsafe {
            ffi::libgap_character_table(
                case as u32,
                slots.as_ptr(),
                slots.len(),
                class_sizes.as_mut_ptr(),
                &mut group_order,
            )
        } == 0
        {
            for slot in &slots {
                unsafe {
                    ffi::libgap_cyc_release(*slot);
                }
            }
            self.slots.borrow_mut().free.extend(slots);
            return Err(error());
        }

        let values: Vec<_> = slots
            .into_iter()
            .map(|slot| Cyclotomic {
                slot,
                slots: Rc::clone(&self.slots),
                _not_send_or_sync: PhantomData,
            })
            .collect();
        let mut values = values.into_iter();
        let rows = (0..row_count)
            .map(|_| values.by_ref().take(column_count).collect())
            .collect();
        Ok(CharacterTable {
            case,
            rows,
            class_sizes,
            group_order,
        })
    }

    pub fn character_tensor_decomposition(
        &self,
        case: CharacterTableCase,
        lhs: usize,
        rhs: usize,
        row_count: usize,
    ) -> Result<Vec<i64>, Error> {
        let mut multiplicities = vec![0; row_count];
        if unsafe {
            ffi::libgap_character_tensor_decomposition(
                case as u32,
                lhs,
                rhs,
                multiplicities.as_mut_ptr(),
                multiplicities.len(),
            )
        } == 0
        {
            return Err(error());
        }
        Ok(multiplicities)
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
    use std::sync::Mutex;

    static LIBGAP_TEST_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn unmodified_gap_handles_rational_cyclotomics() {
        let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
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

    #[test]
    fn unmodified_gap_exposes_character_tables_and_tensor_products() {
        let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
        let context = Context::new().unwrap();
        let table = context.character_table(CharacterTableCase::A5).unwrap();
        assert_eq!(table.group_order, 60);
        assert_eq!(table.class_sizes, [1, 15, 20, 12, 12]);
        assert_eq!(table.rows.len(), 5);
        assert!(table.rows.iter().all(|row| row.len() == 5));
        assert_eq!(
            context
                .character_tensor_decomposition(CharacterTableCase::A5, 1, 1, table.rows.len())
                .unwrap(),
            [1, 1, 1, 1, 1]
        );
    }
}
