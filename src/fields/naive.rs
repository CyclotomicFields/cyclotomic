extern crate num;

use std::cmp::Eq;
use std::ops::{Add, Mul};
use std::vec::Vec;
use self::num::{One, Zero};
use std::ptr::eq;
use crate::fields::FieldElement;

type Z = num::bigint::BigInt;
type Q = num::rational::BigRational;

#[derive(Debug, Clone)]
pub struct Number {
    order: u64     /*   n   */,
    coeffs: Vec<Q> /*   c   */,
    exps: Vec<u64> /*   e   */,
}

impl Number {
    pub fn new(order: u64, coeffs: Vec<Q>, exps: Vec<u64>) -> Number {
        Number {
            order,
            coeffs,
            exps,
        }
    }
}

impl Zero for Number {
    fn zero() -> Self {
        Number::new(1, vec![Q::new(Z::from(0), Z::from(1))], vec![0])
    }

    fn set_zero(&mut self) {
        self.order = 1;
        self.coeffs = vec![Q::new(Z::from(0), Z::from(1))];
        self.exps = vec![]
    }

    fn is_zero(&self) -> bool {
        return self.order == 1
            && self.coeffs == vec![Q::new(Z::from(0), Z::from(1))]
            && self.exps == vec![];
    }
}

impl One for Number {
    fn one() -> Self {
        Number::new(1, vec![Q::new(Z::from(1), Z::from(1))], vec![0])
    }
}

impl Number {
    pub fn increase_order_to(z: &mut Number, new_order: u64) -> () {
        z.exps = z
            .exps
            .clone()
            .into_iter()
            .map(|k| new_order * k / z.order)
            .collect();

        z.order = new_order
    }

    /// TODO: Use Rob's code to actually do some reductions here?
    pub fn match_orders(z1: &mut Number, z2: &mut Number) -> () {
        let new_order = num::integer::lcm(z1.order, z2.order);
        Number::increase_order_to(z1, new_order);
        Number::increase_order_to(z2, new_order);
    }
}

impl PartialEq for Number {
    fn eq(&self, other: &Self) -> bool {
        let mut z1 = self.clone();
        let mut z2 = other.clone();
        Number::match_orders(&mut z1, &mut z2);
        // Since we assume the exponent list is sorted, we can just do a
        // comparison of the coeff and exp vectors, the number is in
        // canonical form.
        z1.coeffs == z2.coeffs && z1.exps == z2.exps
    }

    fn ne(&self, other: &Self) -> bool {
        return !eq(self, other);
    }
}

impl Eq for Number {}

impl Add for Number {
    type Output = Number;

    fn add(self, rhs: Self) -> Self::Output {
        let mut z1 = self.clone();
        let mut z2 = rhs.clone();
        Number::match_orders(&mut z1, &mut z2);

        let mut z1_i = 0;
        let mut z2_i = 0;
        let mut result = Number::zero();

        // We are producing a result cyclotomic in canonical form assuming the
        // inputs are in canonical form. That is, we build the result term by
        // term, assuming that the exponent lists are sorted.
        while z1_i < z1.exps.len() || z2_i < z2.exps.len() {
            while z1_i < z1.exps.len() && z1.exps[z1_i] < z2.exps[z2_i] {
                result.coeffs.push(z1.coeffs[z1_i].clone());
                result.exps.push(z1.exps[z1_i].clone());
                z1_i += 1;
            }

            while z2_i < z2.exps.len() && z2.exps[z2_i] < z1.exps[z1_i] {
                result.coeffs.push(z2.coeffs[z2_i].clone());
                result.exps.push(z2.exps[z2_i].clone());
                z2_i += 1;
            }

            // If the exponents match up, we just add the coefficients - they are for the same term
            if z1.exps[z1_i] == z2.exps[z2_i] {
                result
                    .coeffs
                    .push(z1.coeffs[z1_i].clone() + z2.coeffs[z2_i].clone());
                result.exps.push(z1.exps[z1_i]);
                z1_i += 1;
                z2_i += 1
            }
        }

        result
    }
}

impl Mul for Number {
    type Output = Number;

    /// I could have a guess at how to implement this, but I would inevitably embarrass myself.
    fn mul(self, rhs: Self) -> Self::Output {
        // placeholder?
        self
    }
}

impl FieldElement for Number {
    /// Intentionally bad, makes no effort to avoid copies.
    fn add_mut(&mut self, other: &mut Self) -> () {
        let result = self.clone() + other.clone();
        *self = result;
    }

    /// Also intentionally bad.
    fn mul_mut(&mut self, other: &mut Self) -> () {
        let result = self.clone() * other.clone();
        *self = result;
    }

    // TODO: this is wrong lmao
    fn inv(&self) -> Self {
        self.clone()
    }

    fn inv_mut(&mut self) {
        let result = self.inv();
        *self = result;
    }
}