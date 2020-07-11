extern crate num;
extern crate rustc_hash;

use self::num::{One, Zero};
use crate::fields::{AdditiveGroup, MultiplicativeGroup};
use crate::fields::{CyclotomicFieldElement, FieldElement, Q, Z};
use basis::convert_to_base;
use num::traits::Inv;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::TryInto;
use std::fmt;
use std::hash::BuildHasherDefault;
use std::hash::Hasher;
use std::ops::{AddAssign, Mul, SubAssign};
use std::vec::Vec;
use std::intrinsics::transmute;
use rustc_hash::FxHashMap;
use crate::fields::sparse::basis::try_reduce;

pub mod add;
pub mod basis;
pub mod galois;
pub mod mul;


type ExpCoeffMap = FxHashMap<i64, Q>;

/// Represents a polynomial in the `order`th root of unity.
#[derive(Clone)]
pub struct Number {
    order: i64,
    coeffs: ExpCoeffMap,
}

pub fn print_gap(z: &Number) -> String {
    let mut str_list: Vec<String> = vec![];
    for exp in 0..z.order {
        let zero = Q::zero().clone();
        let coeff = z.coeffs.get(&exp).unwrap_or(&zero);
        if !coeff.is_zero() {
            str_list.push(String::from(
                format!("{} * E({})^{}", coeff, z.order, exp).as_str(),
            ))
        }
    }
    "(".to_string() + &str_list.join(" + ") + ")"
}

impl fmt::Debug for Number {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Number ({})", print_gap(self))
    }
}

impl Number {
    pub fn new(order: i64, coeffs: &ExpCoeffMap) -> Number {
        Number {
            order: order,
            coeffs: coeffs.clone(),
        }
    }

    pub fn increase_order_to(z: &mut Self, new_order: i64) {
        let mut new_coeffs = ExpCoeffMap::default();
        for (exp, coeff) in &z.coeffs {
            new_coeffs.insert(new_order * exp.clone() / z.order.clone(), coeff.clone());
        }
        z.order = new_order.clone();
        z.coeffs = new_coeffs;
    }

    pub fn match_orders(z1: &mut Number, z2: &mut Number) {
        if z1.order == z2.order {
            return;
        }
        let new_order = num::integer::lcm(z1.order, z2.order);
        Number::increase_order_to(z1, new_order);
        Number::increase_order_to(z2, new_order);
    }
}

fn are_coprime(x: i64, y: i64) -> bool {
    num::integer::gcd(x as u64, y as u64) == 1
}

fn phi(n: i64) -> i64 {
    let mut count = 0;
    for k in 1..n {
        if are_coprime(n, k) {
            count += 1;
        }
    }
    count
}

fn get_same_coeff(z: &Number) -> Option<Q> {
    let coeffs = z.coeffs.clone().into_iter().map(|(exp, coeff)| coeff);
    let nonzero_coeffs: HashSet<Q> = coeffs.filter(|q| *q != Q::zero()).collect();

    if nonzero_coeffs.len() == 0 {
        // all coeffs are zero
        Some(Q::zero())
    } else if nonzero_coeffs.len() == 1 {
        Some(nonzero_coeffs.iter().last().unwrap().clone())
    } else {
        None
    }
}

fn math_mod(x: &i64, n: &i64) -> i64 {
    let res = (x % n + n) % n;
    res
}

#[derive(Eq, PartialEq)]
enum Sign {
    Plus,
    Minus,
}

fn add_single(coeffs: &mut ExpCoeffMap, exp: i64, coeff: &Q, sign: Sign) {
    let maybe_existing_coeff = coeffs.get_mut(&exp);
    match maybe_existing_coeff {
        None => {
            if sign == Sign::Plus {
                coeffs.insert(exp, coeff.clone());
            } else {
                coeffs.insert(exp, -coeff.clone());
            }
        }
        Some(existing_coeff) => {
            // TODO: find a way to get rid of coeff.clone() here, it's not needed
            if sign == Sign::Plus {
                existing_coeff.add_assign(coeff.clone());
            } else {
                existing_coeff.sub_assign(coeff.clone());
            }
        }
    }
}

pub fn is_zero(z: &Number) -> bool {
    for (_, coeff) in &z.coeffs {
        if !coeff.is_zero() {
            return false;
        }
    }
    true
}

fn count_powers(n: &i64, n_divisors: &Vec<i64>) -> Vec<(i64, i64)> {
    let mut result = vec![];
    let mut n_factored = n.clone();

    for divisor in n_divisors {
        let mut power: u64 = 0;

        while n_factored % divisor == 0 {
            power += 1;
            n_factored = n_factored / divisor;
        }

        if power != 0 {
            result.push((divisor.clone(), power as i64));
        }
    }

    result
}

impl FieldElement for Number {
    fn eq(&mut self, other: &mut Self) -> bool {
        let mut za = self.clone();
        let mut zb = other.clone();
        Number::match_orders(&mut za, &mut zb);
        let mut z1 = convert_to_base(&za);
        let mut z2 = convert_to_base(&zb);

        // Now that we've matched the orders, z1 and z2 are expressed as
        // elements in the same field so are the same iff each nonzero term is
        // the same.
        fn has_diff(left: &Number, right: &Number) -> bool {
            for (exp_left, coeff_left) in &left.coeffs {
                match right.coeffs.get(&exp_left) {
                    None => {
                        if coeff_left != &Q::zero() {
                            return true;
                        }
                    }
                    Some(coeff_right) => {
                        if coeff_left != coeff_right {
                            return true;
                        }
                    }
                }
            }
            false
        }

        !has_diff(&z1, &z2) && !has_diff(&z2, &z1)
    }
}

impl CyclotomicFieldElement for Number {
    fn e(n: i64, k: i64) -> Self {
        Number::new(
            n,
            &[(k, Q::from_integer(Z::one()))].iter().cloned().collect(),
        )
    }

    fn scalar_mul(&mut self, scalar: &Q) -> &mut Self {
        let mut result = self.clone();
        for (_, coeff) in result.coeffs.iter_mut() {
            *coeff *= scalar.clone();
        }
        *self = result;
        self
    }

    fn zero_order(n: i64) -> Number {
        Number::new(n, &ExpCoeffMap::default())
    }

    fn one_order(n: i64) -> Number {
        let mut coeffs = ExpCoeffMap::default();
        for i in 1..n {
            coeffs.insert(i, -Q::one());
        }
        Number::new(n, &coeffs)
    }
}

pub fn random_rational<G>(g: &mut G) -> Q
where
    G: rand::RngCore,
{
    let p: i64 = g.gen_range(1, 10);
    let q: i64 = g.gen_range(1, 10);
    Q::new(Z::from(p), Z::from(q))
}

pub fn random_cyclotomic<G>(g: &mut G, min_order: i64, max_order: i64) -> Number
where
    G: rand::RngCore,
{
    let order = g.gen_range(min_order, max_order);
    let num_terms: u64 = g.gen_range(1, 5);
    let mut result = Number::zero_order(order.clone());

    for _ in 1..=num_terms {
        let exp: i64 = g.gen_range(1, order);
        let coeff = random_rational(g);
        result.coeffs.insert(exp, coeff);
    }

    result
}

// These are precisely the field axioms.
// TODO: write these tests for the trait, then instantiate somehow? can rust do that?
#[cfg(test)]
mod tests {
    use super::*;

    quickcheck! {
    fn zero_is_add_identity(z: Number) -> bool {
        z.clone().add(&mut Number::zero_order(z.order.clone())).eq(&mut z.clone())
    }
    }

    quickcheck! {
    fn add_is_associative(x: Number, y: Number, z: Number) -> bool {
        (x.clone().add(&mut y.clone())).add(&mut z.clone()).eq(x.clone().add(y.clone().add(&mut z.clone())))
    }
    }

    quickcheck! {
    fn add_is_commutative(x: Number, y: Number) -> bool {
        x.clone().add(&mut y.clone()).eq(y.clone().add(&mut x.clone()))
    }
    }

    quickcheck! {
    fn one_is_mul_identity(z: Number) -> bool {
        let mut same = z.clone().mul(&mut Number::one_order(z.order.clone())).clone();
        same.eq(&mut z.clone())
    }
    }

    quickcheck! {
    fn add_has_inverses(z: Number) -> bool {
        z.clone().add(z.clone().add_invert()).eq(&mut Number::zero_order(z.order))
    }
    }

    quickcheck! {
    fn zero_kills_all(z: Number) -> bool {
        Number::zero_order(z.order.clone()).mul(&mut z.clone()).eq(&mut Number::zero_order(z.order))
    }
    }

    quickcheck! {
    fn mul_is_commutative(x: Number, y: Number) -> bool {
        x.clone().mul(&mut y.clone()).eq(y.clone().mul(&mut x.clone()))
    }
    }

    quickcheck! {
    fn mul_is_associative(x: Number, y: Number, z: Number) -> bool {
        (x.clone().mul(&mut y.clone())).mul(&mut z.clone()).eq(x.clone().mul(y.clone().mul(&mut z.clone())))
    }
    }

    quickcheck! {
    fn mul_has_inverses(arb_z: Number) -> bool {
        let z = convert_to_base(&arb_z);
        if is_zero(&z) {
            return true;
        }
        let mut prod = z.clone().mul_invert().mul(&mut z.clone()).clone();
        prod.eq(&mut Number::one_order(z.order))
    }
    }

    quickcheck! {
    fn mul_distributes_over_add(x: Number, y: Number, z: Number) -> bool {
        x.clone().mul(y.clone().add(&mut z.clone())).eq(x.clone().mul(&mut y.clone()).add(x.clone().mul(&mut z.clone())))
    }
    }

    impl Arbitrary for Number {
        fn arbitrary<G>(g: &mut G) -> Self
        where
            G: Gen,
        {
            random_cyclotomic(g, 100, 100)
        }
    }
}
