extern crate num;

use self::num::{One, Zero};
use crate::fields::{CyclotomicFieldElement, Exponent, FieldElement, Q, Z};
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use std::cmp::Eq;
use std::collections::HashMap;
use std::fmt;
use std::ops::{Add, Mul};
use std::ptr::eq;
use std::vec::Vec;

/// Represents a polynomial in the `order`th root of unity.
///
/// Simplest possible choice of data structure - a hash map. We assume only that
/// each term has an exponent less than the order. So the only real term is the
/// one for exponent zero, for example.
#[derive(Clone)]
pub struct Number {
    order: Exponent,
    coeffs: HashMap<Exponent, Q>,
}

impl fmt::Debug for Number {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut str_list: Vec<String> = vec![];
        for exp in 0..self.order {
            let z = Q::zero().clone();
            let coeff = self.coeffs.get(&exp).unwrap_or(&z);
            if *coeff != z {
                str_list.push(String::from(
                    format!("{} * E({})^{}", coeff, self.order, exp).as_str(),
                ))
            }
        }

        write!(f, "Number ({})", str_list.join(" + "))
    }
}

impl Number {
    pub fn new(order: Exponent, coeffs: HashMap<Exponent, Q>) -> Number {
        Number { order, coeffs }
    }
}

/// Represents zero iff all coeffs are zero. Note that the encoding of
/// zero is far from unique.
impl Zero for Number {
    fn zero() -> Self {
        Number::new(
            1,
            [(0, Q::from_integer(Z::zero()))].iter().cloned().collect(),
        )
    }

    fn set_zero(&mut self) {
        self.order = 1;
        self.coeffs = [(0, Q::from_integer(Z::zero()))].iter().cloned().collect();
    }

    fn is_zero(&self) -> bool {
        self.coeffs.values().all(|coeff| coeff.is_zero())
    }
}

impl One for Number {
    fn one() -> Self {
        Number::new(
            1,
            [(0, Q::from_integer(Z::one()))].iter().cloned().collect(),
        )
    }
}

impl Number {
    pub fn increase_order_to(z: &mut Self, new_order: u64) -> () {
        let mut new_coeffs = HashMap::new();
        for (exp, coeff) in z.coeffs.clone() {
            new_coeffs.insert(new_order * exp / z.order, coeff);
        }
        z.order = new_order;
        z.coeffs = new_coeffs;
    }

    /// TODO: Use Rob's code to actually do some reductions here?
    pub fn match_orders(z1: &mut Number, z2: &mut Number) -> () {
        let new_order = num::integer::lcm(z1.order, z2.order);
        Number::increase_order_to(z1, new_order);
        Number::increase_order_to(z2, new_order);
        assert_eq!(z1.order, z2.order);
    }
}

impl PartialEq for Number {
    fn eq(&self, other: &Self) -> bool {
        let mut z1 = self.clone();
        let mut z2 = other.clone();
        Number::match_orders(&mut z1, &mut z2);

        // Now that we've matched the orders, z1 and z2 are expressed as
        // elements in the same field, so are the same iff each term is the
        // same.
        z1.coeffs == z2.coeffs
    }

    fn ne(&self, other: &Self) -> bool {
        return !eq(self, other);
    }
}

impl Eq for Number {}

impl Add for Number {
    type Output = Number;

    /// Simplest possible - term wise addition using hashing.
    ///
    /// Purposely written so it is obviously symmetric in the parameters, thus
    /// commutative by inspection. Of course, there are tests for that.
    fn add(self, rhs: Self) -> Self::Output {
        let mut z1 = self.clone();
        let mut z2 = rhs.clone();
        Self::match_orders(&mut z1, &mut z2);

        // We will never need to reduce here, you can't add low powers of
        // $\zeta_n$ and get higher powers. Higher powers do not exist.

        let mut coeffs: HashMap<u64, Q> = HashMap::new();
        for (exp, coeff) in z1.coeffs.into_iter().chain(z2.coeffs) {
            match coeffs.get(&exp) {
                Some(existing_coeff) => coeffs.insert(exp, coeff + existing_coeff),
                None => coeffs.insert(exp, coeff),
            };
        }

        let result = Number::new(z1.order, coeffs);
        result
    }
}

impl Mul for Number {
    type Output = Number;

    /// Multiplies term by term, not bothering to do anything interesting.
    fn mul(self, rhs: Self) -> Self::Output {
        let mut z1 = self.clone();
        let mut z2 = rhs.clone();
        Self::match_orders(&mut z1, &mut z2);

        let mut result = Self::zero();

        // This order is almost certainly not optimal. But you know, whatever.
        result.order = z1.order;
        for (exp1, coeff1) in z1.coeffs.clone() {
            for (exp2, coeff2) in z2.coeffs.clone() {
                // The mod order part is the only reduction we do.
                // $\zeta_n^n = 1$ is still trivial enough to include in a naive
                // implementation, I think.
                let new_exp = (exp1 + exp2) % z1.order;
                let new_coeff = coeff1.clone() * coeff2.clone();
                result.coeffs.insert(new_exp, new_coeff);
            }
        }
        result
    }
}

impl FieldElement for Number {
    /// Intentionally bad, makes no effort to avoid copies.
    fn add_mut(&mut self, other: &mut Self) -> Self {
        self.clone() + other.clone()
    }

    /// Also intentionally bad.
    fn mul_mut(&mut self, other: &mut Self) -> Self {
        self.clone() * other.clone()
    }

    /// TODO: this is wrong lmao, need polynomial division? or is there a faster way?
    fn inv(&self) -> Self {
        self.clone()
    }

    /// TODO: Also wrong. Also lmao.
    fn inv_mut(&mut self) -> Self {
        self.inv()
    }
}

impl CyclotomicFieldElement for Number {
    fn e(n: u64, k: u64) -> Self {
        Number::new(
            n,
            [(k, Q::from_integer(Z::one()))].iter().cloned().collect(),
        )
    }

    fn scalar_mul(&self, scalar: Q) -> Self {
        let mut result = self.clone();
        for (_, coeff) in result.coeffs.iter_mut() {
            *coeff *= scalar.clone();
        }
        result
    }
}

// Just a hack to get arbitrary coefficients
#[derive(Clone)]
struct QArb(Q);

impl Arbitrary for QArb {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        let p: i64 = g.gen_range(-1024, 1024);
        let q: i64 = g.gen_range(1, 1024);
        QArb(Q::new(Z::from(p), Z::from(q)))
    }
}

impl Arbitrary for Number {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        let order: u64 = g.gen_range(1, 1000);
        let num_terms: u64 = g.gen_range(1, 10);
        let mut result = Self::zero();
        result.order = order;

        for _ in 1..=num_terms {
            let exp: u64 = g.gen_range(0, order);
            let QArb(coeff) = QArb::arbitrary(g);
            result.coeffs.insert(exp, coeff);
        }

        result
    }
}

// These are precisely the field axioms.
// TODO: add the other axioms
// TODO: add multiplicative inverses to this
// TODO: write these tests for the trait, then instantiate somehow? can rust do that?
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_neq_one() {
        assert_ne!(Number::zero(), Number::one());
    }

    quickcheck! {
    fn zero_is_add_identity(z: Number) -> bool {
        z.clone() + Number::zero() == z.clone()
    }
    }

    quickcheck! {
    fn one_is_mul_identity(z: Number) -> bool {
        z.clone() * Number::one() == z.clone()
    }
    }

    quickcheck! {
    fn add_is_associative(x: Number, y: Number, z: Number) -> bool {
        (x.clone() + y.clone()) + z.clone() == x.clone() + (y.clone() + z.clone())
    }
    }

    quickcheck! {
    fn add_is_commutative(x: Number, y: Number) -> bool {
        x.clone() + y.clone() == y.clone() + x.clone()
    }
    }
}
