//! Packed GAP-style cyclotomics with inline integers and exact rational
//! promotion.

use crate::literal::{Coefficient, KernelContext, KernelCyclotomic};
use crate::Error;
use rug::{Integer, Rational};
use std::fmt;

/// One-word exact coefficient, analogous to GAP's tagged `Obj`.
///
/// The low bit distinguishes an immediate signed integer from an aligned
/// pointer to a boxed GMP integer or rational. Like GAP immediate integers,
/// the inline range is one bit narrower than the machine word; larger integers
/// promote without entering rational arithmetic.
pub struct TaggedRational(usize);

impl TaggedRational {
    const TAG_SMALL: usize = 1;
    const TAG_INTEGER: usize = 2;
    const TAG_MASK: usize = 3;

    fn small_value(&self) -> Option<i64> {
        if self.0 & Self::TAG_SMALL == 0 {
            None
        } else {
            Some(((self.0 as isize) >> 1) as i64)
        }
    }

    fn encode_small(value: i64) -> Option<usize> {
        let value = isize::try_from(value).ok()?;
        if value < isize::MIN / 2 || value > isize::MAX / 2 {
            return None;
        }
        Some(((value as usize) << 1) | Self::TAG_SMALL)
    }

    fn big_ref(&self) -> &Rational {
        debug_assert_eq!(self.0 & Self::TAG_MASK, 0);
        // SAFETY: every untagged word is created with Box::into_raw below,
        // remains owned by this value, and is freed exactly once in Drop.
        unsafe { &*(self.0 as *const Rational) }
    }

    fn big_mut(&mut self) -> &mut Rational {
        debug_assert_eq!(self.0 & Self::TAG_MASK, 0);
        // SAFETY: coefficients deep-clone boxed rationals, so this value is
        // the unique owner of the pointed-to allocation.
        unsafe { &mut *(self.0 as *mut Rational) }
    }

    fn integer_ref(&self) -> &Integer {
        debug_assert_eq!(self.0 & Self::TAG_MASK, Self::TAG_INTEGER);
        let pointer = self.0 & !Self::TAG_MASK;
        // SAFETY: tagged integer pointers come exclusively from from_integer.
        unsafe { &*(pointer as *const Integer) }
    }

    fn integer_mut(&mut self) -> &mut Integer {
        debug_assert_eq!(self.0 & Self::TAG_MASK, Self::TAG_INTEGER);
        let pointer = self.0 & !Self::TAG_MASK;
        // SAFETY: Clone creates a distinct allocation, so ownership is unique.
        unsafe { &mut *(pointer as *mut Integer) }
    }

    fn is_integer_pointer(&self) -> bool {
        self.0 & Self::TAG_MASK == Self::TAG_INTEGER
    }

    fn is_integer(&self) -> bool {
        self.is_small() || self.is_integer_pointer()
    }

    pub fn from_fraction(numerator: i64, denominator: i64) -> Self {
        assert_ne!(denominator, 0, "rational denominator must be nonzero");
        if denominator == 1 {
            return Self::from(numerator);
        }
        Self::from_big(Rational::from((numerator, denominator)))
    }

    fn from_big(value: Rational) -> Self {
        if value.is_integer() {
            return Self::from_integer(value.numer().clone());
        }
        let pointer = Box::into_raw(Box::new(value)) as usize;
        debug_assert_ne!(pointer, 0);
        debug_assert_eq!(pointer & Self::TAG_SMALL, 0);
        Self(pointer)
    }

    fn from_integer(value: Integer) -> Self {
        if let Some(value) = value.to_i64() {
            if let Some(encoded) = Self::encode_small(value) {
                return Self(encoded);
            }
        }
        let pointer = Box::into_raw(Box::new(value)) as usize;
        debug_assert_ne!(pointer, 0);
        debug_assert_eq!(pointer & Self::TAG_MASK, 0);
        Self(pointer | Self::TAG_INTEGER)
    }

    fn to_integer(&self) -> Integer {
        match self.small_value() {
            Some(value) => Integer::from(value),
            None if self.is_integer_pointer() => self.integer_ref().clone(),
            None => unreachable!("rational coefficient is not an integer"),
        }
    }

    fn to_big(&self) -> Rational {
        match self.small_value() {
            Some(value) => Rational::from(value),
            None if self.is_integer_pointer() => Rational::from(self.integer_ref().clone()),
            None => self.big_ref().clone(),
        }
    }

    pub fn is_small(&self) -> bool {
        self.0 & Self::TAG_SMALL != 0
    }

    pub fn numerator(&self) -> Integer {
        match self.small_value() {
            Some(value) => Integer::from(value),
            None if self.is_integer_pointer() => self.integer_ref().clone(),
            None => self.big_ref().numer().clone(),
        }
    }

    pub fn denominator(&self) -> Integer {
        match self.small_value() {
            Some(_) => Integer::from(1),
            None if self.is_integer_pointer() => Integer::from(1),
            None => self.big_ref().denom().clone(),
        }
    }

    pub fn as_rational(&self) -> Rational {
        self.to_big()
    }
}

impl From<i64> for TaggedRational {
    fn from(value: i64) -> Self {
        match Self::encode_small(value) {
            Some(encoded) => Self(encoded),
            None => Self::from_integer(Integer::from(value)),
        }
    }
}

impl Clone for TaggedRational {
    fn clone(&self) -> Self {
        match self.small_value() {
            Some(_) => Self(self.0),
            None if self.is_integer_pointer() => Self::from_integer(self.integer_ref().clone()),
            None => Self::from_big(self.big_ref().clone()),
        }
    }
}

impl Drop for TaggedRational {
    fn drop(&mut self) {
        if self.is_integer_pointer() {
            let pointer = self.0 & !Self::TAG_MASK;
            // SAFETY: this is the unique pointer allocated in from_integer.
            unsafe {
                drop(Box::from_raw(pointer as *mut Integer));
            }
        } else if self.0 & Self::TAG_SMALL == 0 {
            // SAFETY: this is the unique pointer created by Box::into_raw in
            // from_big. Clone allocates a distinct box.
            unsafe {
                drop(Box::from_raw(self.0 as *mut Rational));
            }
        }
    }
}

impl PartialEq for TaggedRational {
    fn eq(&self, other: &Self) -> bool {
        match (self.small_value(), other.small_value()) {
            (Some(left), Some(right)) => left == right,
            (None, None) if self.is_integer_pointer() && other.is_integer_pointer() => {
                self.integer_ref() == other.integer_ref()
            }
            _ => self.to_big() == other.to_big(),
        }
    }
}

impl Eq for TaggedRational {}

impl fmt::Debug for TaggedRational {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.small_value() {
            Some(value) => fmt::Debug::fmt(&value, formatter),
            None if self.is_integer_pointer() => fmt::Debug::fmt(self.integer_ref(), formatter),
            None => fmt::Debug::fmt(self.big_ref(), formatter),
        }
    }
}

impl Coefficient for TaggedRational {
    fn zero() -> Self {
        Self::from(0)
    }

    fn one() -> Self {
        Self::from(1)
    }

    fn is_zero(&self) -> bool {
        match self.small_value() {
            Some(value) => value == 0,
            None if self.is_integer_pointer() => self.integer_ref() == &0,
            None => self.big_ref().numer() == &0,
        }
    }

    fn is_one(&self) -> bool {
        self.small_value() == Some(1)
    }

    fn is_minus_one(&self) -> bool {
        self.small_value() == Some(-1)
    }

    fn add(&self, rhs: &Self) -> Result<Self, Error> {
        if let (Some(left), Some(right)) = (self.small_value(), rhs.small_value()) {
            if let Some(value) = left.checked_add(right) {
                return Ok(Self::from(value));
            }
        }
        if self.is_integer() && rhs.is_integer() {
            return Ok(Self::from_integer(self.to_integer() + rhs.to_integer()));
        }
        let mut value = self.to_big();
        match rhs.small_value() {
            Some(rhs) => value += rhs,
            None if rhs.is_integer_pointer() => value += rhs.integer_ref(),
            None => value += rhs.big_ref(),
        }
        Ok(Self::from_big(value))
    }

    fn sub(&self, rhs: &Self) -> Result<Self, Error> {
        if let (Some(left), Some(right)) = (self.small_value(), rhs.small_value()) {
            if let Some(value) = left.checked_sub(right) {
                return Ok(Self::from(value));
            }
        }
        if self.is_integer() && rhs.is_integer() {
            return Ok(Self::from_integer(self.to_integer() - rhs.to_integer()));
        }
        let mut value = self.to_big();
        match rhs.small_value() {
            Some(rhs) => value -= rhs,
            None if rhs.is_integer_pointer() => value -= rhs.integer_ref(),
            None => value -= rhs.big_ref(),
        }
        Ok(Self::from_big(value))
    }

    fn add_assign(&mut self, rhs: &Self) -> Result<(), Error> {
        if let (Some(left), Some(right)) = (self.small_value(), rhs.small_value()) {
            if let Some(value) = left.checked_add(right) {
                if let Some(encoded) = Self::encode_small(value) {
                    self.0 = encoded;
                    return Ok(());
                }
            }
            *self = Self::from_integer(Integer::from(left) + right);
            return Ok(());
        }
        if self.is_integer_pointer() && rhs.is_integer() {
            match rhs.small_value() {
                Some(rhs) => *self.integer_mut() += rhs,
                None => *self.integer_mut() += rhs.integer_ref(),
            }
            return Ok(());
        }
        if !self.is_integer() {
            match rhs.small_value() {
                Some(rhs) => *self.big_mut() += rhs,
                None if rhs.is_integer_pointer() => *self.big_mut() += rhs.integer_ref(),
                None => *self.big_mut() += rhs.big_ref(),
            }
            return Ok(());
        }
        *self = if rhs.is_integer() {
            Self::from_integer(self.to_integer() + rhs.to_integer())
        } else {
            Self::from_big(self.to_big() + rhs.big_ref())
        };
        Ok(())
    }

    fn sub_assign(&mut self, rhs: &Self) -> Result<(), Error> {
        if let (Some(left), Some(right)) = (self.small_value(), rhs.small_value()) {
            if let Some(value) = left.checked_sub(right) {
                if let Some(encoded) = Self::encode_small(value) {
                    self.0 = encoded;
                    return Ok(());
                }
            }
            *self = Self::from_integer(Integer::from(left) - right);
            return Ok(());
        }
        if self.is_integer_pointer() && rhs.is_integer() {
            match rhs.small_value() {
                Some(rhs) => *self.integer_mut() -= rhs,
                None => *self.integer_mut() -= rhs.integer_ref(),
            }
            return Ok(());
        }
        if !self.is_integer() {
            match rhs.small_value() {
                Some(rhs) => *self.big_mut() -= rhs,
                None if rhs.is_integer_pointer() => *self.big_mut() -= rhs.integer_ref(),
                None => *self.big_mut() -= rhs.big_ref(),
            }
            return Ok(());
        }
        *self = if rhs.is_integer() {
            Self::from_integer(self.to_integer() - rhs.to_integer())
        } else {
            Self::from_big(self.to_big() - rhs.big_ref())
        };
        Ok(())
    }

    fn mul(&self, rhs: &Self) -> Result<Self, Error> {
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::zero());
        }
        if self.is_one() {
            return Ok(rhs.clone());
        }
        if rhs.is_one() {
            return Ok(self.clone());
        }
        if self.is_minus_one() {
            return rhs.neg();
        }
        if rhs.is_minus_one() {
            return self.neg();
        }
        if let (Some(left), Some(right)) = (self.small_value(), rhs.small_value()) {
            if let Some(value) = left.checked_mul(right) {
                return Ok(Self::from(value));
            }
        }
        if self.is_integer() && rhs.is_integer() {
            return Ok(Self::from_integer(self.to_integer() * rhs.to_integer()));
        }
        let mut value = self.to_big();
        match rhs.small_value() {
            Some(rhs) => value *= rhs,
            None if rhs.is_integer_pointer() => value *= rhs.integer_ref(),
            None => value *= rhs.big_ref(),
        }
        Ok(Self::from_big(value))
    }

    fn neg(&self) -> Result<Self, Error> {
        match self.small_value() {
            Some(value) => match value.checked_neg() {
                Some(value) => Ok(Self::from(value)),
                None => Ok(Self::from_integer(-self.to_integer())),
            },
            None if self.is_integer_pointer() => {
                Ok(Self::from_integer(-self.integer_ref().clone()))
            }
            None => Ok(Self::from_big(-self.to_big())),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Cyclotomic {
    /// All coefficients fit the immediate integer path.
    Small(KernelCyclotomic<i64>),
    /// At least one operation required exact rational/GMP storage.
    Exact(KernelCyclotomic<TaggedRational>),
}

pub struct Context {
    small_kernel: KernelContext<i64>,
    exact_kernel: KernelContext<TaggedRational>,
}

impl Default for Context {
    fn default() -> Self {
        Self {
            small_kernel: KernelContext::default(),
            exact_kernel: KernelContext::default(),
        }
    }
}

impl Context {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_integer_terms(
        &mut self,
        order: u32,
        terms: &[(u32, i64)],
    ) -> Result<Cyclotomic, Error> {
        if let Ok(value) = self.small_kernel.from_terms(order, terms) {
            return Ok(Cyclotomic::Small(value));
        }
        let terms: Vec<_> = terms
            .iter()
            .map(|&(exponent, coefficient)| (exponent, TaggedRational::from(coefficient)))
            .collect();
        self.exact_kernel
            .from_terms(order, &terms)
            .map(Cyclotomic::Exact)
    }

    pub fn from_terms(
        &mut self,
        order: u32,
        terms: &[(u32, (i64, i64))],
    ) -> Result<Cyclotomic, Error> {
        if terms.iter().any(|(_, (_, denominator))| *denominator == 0) {
            return Err(Error("rational denominator must be nonzero".into()));
        }
        if terms.iter().all(|(_, (_, denominator))| *denominator == 1) {
            let integer_terms: Vec<_> = terms
                .iter()
                .map(|&(exponent, (numerator, _))| (exponent, numerator))
                .collect();
            return self.from_integer_terms(order, &integer_terms);
        }
        let mut common_denominator = 1_i128;
        let mut valid_common_denominator = true;
        for (_, (_, denominator)) in terms {
            let denominator = i128::from(*denominator).abs();
            let mut left = common_denominator;
            let mut right = denominator;
            while right != 0 {
                (left, right) = (right, left % right);
            }
            match common_denominator
                .checked_div(left)
                .and_then(|value| value.checked_mul(denominator))
            {
                Some(value) if value <= i128::from(i64::MAX) => common_denominator = value,
                _ => {
                    valid_common_denominator = false;
                    break;
                }
            }
        }
        if valid_common_denominator {
            let integer_terms: Option<Vec<_>> = terms
                .iter()
                .map(|&(exponent, (numerator, denominator))| {
                    let sign = if denominator < 0 { -1_i128 } else { 1_i128 };
                    let multiplier = common_denominator / i128::from(denominator).abs();
                    let coefficient = i128::from(numerator)
                        .checked_mul(sign)?
                        .checked_mul(multiplier)?;
                    Some((exponent, i64::try_from(coefficient).ok()?))
                })
                .collect();
            if let Some(integer_terms) = integer_terms {
                let integer = self.from_integer_terms(order, &integer_terms)?;
                return self.scale_fraction(
                    &integer,
                    1,
                    i64::try_from(common_denominator).unwrap(),
                );
            }
        }
        let terms: Vec<_> = terms
            .iter()
            .map(|&(exponent, (numerator, denominator))| {
                (
                    exponent,
                    TaggedRational::from_fraction(numerator, denominator),
                )
            })
            .collect();
        self.exact_kernel
            .from_terms(order, &terms)
            .map(Cyclotomic::Exact)
    }

    pub fn root(&mut self, order: u32, exponent: u32) -> Result<Cyclotomic, Error> {
        self.small_kernel
            .root(order, exponent)
            .map(Cyclotomic::Small)
    }

    pub fn add(&mut self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        if let (Cyclotomic::Small(lhs), Cyclotomic::Small(rhs)) = (lhs, rhs) {
            if let Ok(value) = self.small_kernel.add(lhs, rhs) {
                return Ok(Cyclotomic::Small(value));
            }
        }
        let lhs = Self::promote(lhs);
        let rhs = Self::promote(rhs);
        self.exact_kernel.add(&lhs, &rhs).map(Cyclotomic::Exact)
    }

    pub fn add_assign(&mut self, lhs: &mut Cyclotomic, rhs: &Cyclotomic) -> Result<(), Error> {
        if let (Cyclotomic::Small(left), Cyclotomic::Small(right)) = (&mut *lhs, rhs) {
            if self.small_kernel.add_assign(left, right).is_ok() {
                return Ok(());
            }
        }
        let left = Self::promote(lhs);
        let right = Self::promote(rhs);
        *lhs = Cyclotomic::Exact(self.exact_kernel.add(&left, &right)?);
        Ok(())
    }

    pub fn mul(&mut self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        if let (Cyclotomic::Small(lhs), Cyclotomic::Small(rhs)) = (lhs, rhs) {
            if let Ok(value) = self.small_kernel.mul(lhs, rhs) {
                return Ok(Cyclotomic::Small(value));
            }
        }
        let lhs = Self::promote(lhs);
        let rhs = Self::promote(rhs);
        self.exact_kernel.mul(&lhs, &rhs).map(Cyclotomic::Exact)
    }

    pub fn mul_add_assign(
        &mut self,
        accumulator: &mut Cyclotomic,
        lhs: &Cyclotomic,
        rhs: &Cyclotomic,
    ) -> Result<(), Error> {
        if let (Cyclotomic::Small(accumulator), Cyclotomic::Small(left), Cyclotomic::Small(right)) =
            (&mut *accumulator, lhs, rhs)
        {
            if self
                .small_kernel
                .mul_add_assign(accumulator, left, right)
                .is_ok()
            {
                return Ok(());
            }
        }
        let mut exact_accumulator = Self::promote(accumulator);
        let left = Self::promote(lhs);
        let right = Self::promote(rhs);
        self.exact_kernel
            .mul_add_assign(&mut exact_accumulator, &left, &right)?;
        *accumulator = Cyclotomic::Exact(exact_accumulator);
        Ok(())
    }

    pub fn sum_products(
        &mut self,
        products: &[(&Cyclotomic, &Cyclotomic)],
    ) -> Result<Cyclotomic, Error> {
        if products.iter().all(|(lhs, rhs)| {
            matches!(lhs, Cyclotomic::Small(_)) && matches!(rhs, Cyclotomic::Small(_))
        }) {
            let small_products: Vec<_> = products
                .iter()
                .map(|(lhs, rhs)| match (lhs, rhs) {
                    (Cyclotomic::Small(lhs), Cyclotomic::Small(rhs)) => (lhs, rhs),
                    _ => unreachable!(),
                })
                .collect();
            if let Ok(result) = self.small_kernel.sum_products(&small_products) {
                return Ok(Cyclotomic::Small(result));
            }
        }

        let exact_products: Vec<_> = products
            .iter()
            .map(|(lhs, rhs)| (Self::promote(lhs), Self::promote(rhs)))
            .collect();
        let exact_product_refs: Vec<_> =
            exact_products.iter().map(|(lhs, rhs)| (lhs, rhs)).collect();
        self.exact_kernel
            .sum_products(&exact_product_refs)
            .map(Cyclotomic::Exact)
    }

    pub fn sum_products_integer_quotient(
        &mut self,
        products: &[(&Cyclotomic, &Cyclotomic)],
        divisor: i64,
    ) -> Result<i64, Error> {
        if products.iter().all(|(lhs, rhs)| {
            matches!(lhs, Cyclotomic::Small(_)) && matches!(rhs, Cyclotomic::Small(_))
        }) {
            let small_products: Vec<_> = products
                .iter()
                .map(|(lhs, rhs)| match (lhs, rhs) {
                    (Cyclotomic::Small(lhs), Cyclotomic::Small(rhs)) => (lhs, rhs),
                    _ => unreachable!(),
                })
                .collect();
            if let Ok(result) = self
                .small_kernel
                .sum_products_integer_quotient(&small_products, divisor)
            {
                return Ok(result);
            }
            if let Some((order, coefficients)) =
                self.small_kernel.sum_products_crt_basis(&small_products)?
            {
                let terms: Vec<_> = coefficients
                    .into_iter()
                    .enumerate()
                    .filter(|(_, coefficient)| coefficient != &0)
                    .map(|(exponent, coefficient)| {
                        (exponent as u32, TaggedRational::from_integer(coefficient))
                    })
                    .collect();
                let sum = self.exact_kernel.from_basis_terms(order, &terms)?;
                let sum = Cyclotomic::Exact(sum);
                let quotient = self.scale_fraction(&sum, 1, divisor)?;
                return integer_i64(&quotient);
            }
        }

        let sum = self.sum_products(products)?;
        let quotient = self.scale_fraction(&sum, 1, divisor)?;
        integer_i64(&quotient)
    }

    pub fn sum_product_rows_integer_quotients(
        &mut self,
        products: &[Cyclotomic],
        weights: &[Vec<Cyclotomic>],
        divisor: i64,
    ) -> Result<Vec<i64>, Error> {
        if products
            .iter()
            .all(|value| matches!(value, Cyclotomic::Small(_)))
            && weights
                .iter()
                .flatten()
                .all(|value| matches!(value, Cyclotomic::Small(_)))
        {
            let small_products: Vec<_> = products
                .iter()
                .map(|value| match value {
                    Cyclotomic::Small(value) => value,
                    _ => unreachable!(),
                })
                .collect();
            let small_weights: Vec<_> = weights
                .iter()
                .flatten()
                .map(|value| match value {
                    Cyclotomic::Small(value) => value,
                    _ => unreachable!(),
                })
                .collect();
            if let Ok(results) = self.small_kernel.sum_product_rows_integer_quotients(
                &small_products,
                &small_weights,
                divisor,
            ) {
                return Ok(results);
            }
        }

        weights
            .iter()
            .map(|row| {
                let terms: Vec<_> = products.iter().zip(row).collect();
                self.sum_products_integer_quotient(&terms, divisor)
            })
            .collect()
    }

    pub fn conjugate(&mut self, value: &Cyclotomic) -> Result<Cyclotomic, Error> {
        match value {
            Cyclotomic::Small(value) => self.small_kernel.conjugate(value).map(Cyclotomic::Small),
            Cyclotomic::Exact(value) => self.exact_kernel.conjugate(value).map(Cyclotomic::Exact),
        }
    }

    pub fn scale_integer(&mut self, value: &Cyclotomic, scalar: i64) -> Result<Cyclotomic, Error> {
        if let Cyclotomic::Small(value) = value {
            if let Ok(result) = self.small_kernel.scale(value, &scalar) {
                return Ok(Cyclotomic::Small(result));
            }
        }
        let value = Self::promote(value);
        self.exact_kernel
            .scale(&value, &TaggedRational::from(scalar))
            .map(Cyclotomic::Exact)
    }

    pub fn scale_fraction(
        &mut self,
        value: &Cyclotomic,
        numerator: i64,
        denominator: i64,
    ) -> Result<Cyclotomic, Error> {
        if denominator == 1 {
            return self.scale_integer(value, numerator);
        }
        let value = Self::promote(value);
        self.exact_kernel
            .scale(
                &value,
                &TaggedRational::from_fraction(numerator, denominator),
            )
            .map(Cyclotomic::Exact)
    }

    fn promote(value: &Cyclotomic) -> KernelCyclotomic<TaggedRational> {
        match value {
            Cyclotomic::Small(value) => {
                value.map_coefficients(|coefficient| TaggedRational::from(*coefficient))
            }
            Cyclotomic::Exact(value) => value.clone(),
        }
    }
}

fn integer_i64(value: &Cyclotomic) -> Result<i64, Error> {
    if value.order() != 1 {
        return Err(Error("cyclotomic is not rational".into()));
    }
    let terms = value.rational_terms();
    match terms.as_slice() {
        [] => Ok(0),
        [(0, coefficient)] if coefficient.denom() == &1 => coefficient
            .numer()
            .to_i64()
            .ok_or_else(|| Error("integer does not fit i64".into())),
        _ => Err(Error("cyclotomic is not an integer".into())),
    }
}

impl Cyclotomic {
    pub fn order(&self) -> u32 {
        match self {
            Self::Small(value) => value.order(),
            Self::Exact(value) => value.order(),
        }
    }

    pub fn terms(&self) -> Vec<(u32, TaggedRational)> {
        match self {
            Self::Small(value) => value
                .terms()
                .into_iter()
                .map(|(exponent, coefficient)| (exponent, TaggedRational::from(coefficient)))
                .collect(),
            Self::Exact(value) => value.terms(),
        }
    }

    pub fn rational_terms(&self) -> Vec<(u32, Rational)> {
        match self {
            Self::Small(value) => value
                .terms()
                .into_iter()
                .map(|(exponent, coefficient)| (exponent, Rational::from(coefficient)))
                .collect(),
            Self::Exact(value) => value
                .terms()
                .into_iter()
                .map(|(exponent, coefficient)| (exponent, coefficient.as_rational()))
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coefficients_stay_small_and_promote_exactly() {
        assert_eq!(
            std::mem::size_of::<TaggedRational>(),
            std::mem::size_of::<usize>()
        );
        let small = TaggedRational::from(6)
            .mul(&TaggedRational::from(7))
            .unwrap();
        assert!(small.is_small());
        assert_eq!(small.as_rational(), 42);

        let fraction = TaggedRational::from_fraction(1, 3)
            .add(&TaggedRational::from_fraction(1, 6))
            .unwrap();
        assert!(!fraction.is_small());
        assert_eq!(fraction.as_rational(), Rational::from((1, 2)));

        let overflow = TaggedRational::from(i64::MAX)
            .add(&TaggedRational::from(1))
            .unwrap();
        assert!(!overflow.is_small());
        assert!(overflow.is_integer_pointer());
        assert_eq!(overflow.numerator(), Integer::from(i64::MAX) + 1);
        assert_eq!(overflow.denominator(), 1);
        assert!(!fraction.is_integer_pointer());
    }

    #[test]
    fn integer_quotient_promotes_on_intermediate_overflow() {
        let mut context = Context::new();
        let maximum = context.from_terms(1, &[(0, (i64::MAX, 1))]).unwrap();
        let two = context.from_terms(1, &[(0, (2, 1))]).unwrap();

        assert_eq!(
            context
                .sum_products_integer_quotient(&[(&maximum, &two)], 2)
                .unwrap(),
            i64::MAX
        );
    }

    #[test]
    fn multi_output_integer_quotient_preserves_overflow_fallback() {
        let mut context = Context::new();
        let maximum = context.from_terms(1, &[(0, (i64::MAX, 1))]).unwrap();
        let two = context.from_terms(1, &[(0, (2, 1))]).unwrap();

        assert_eq!(
            context
                .sum_product_rows_integer_quotients(
                    std::slice::from_ref(&maximum),
                    &[vec![two]],
                    2,
                )
                .unwrap(),
            vec![i64::MAX]
        );
    }
}
