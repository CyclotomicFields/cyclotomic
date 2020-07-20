extern crate num;
extern crate rustc_hash;

use self::num::{One, Zero};
use crate::fields::sparse_vec::basis::try_reduce;
use crate::fields::{AdditiveGroupElement, MultiplicativeGroupElement};
use crate::fields::{CyclotomicFieldElement, FieldElement, Q, Z};
use basis::convert_to_base;
use num::traits::Inv;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::TryInto;
use std::fmt;
use std::hash::BuildHasherDefault;
use std::hash::Hasher;
use std::intrinsics::transmute;
use std::ops::{AddAssign, Mul, SubAssign};
use std::vec::Vec;

#[macro_use]
use crate::fields::*;

pub mod add;
pub mod basis;
pub mod galois;
pub mod mul;

/// Represents a polynomial in the `order`th root of unity.
#[derive(Clone)]
pub struct Number {
    order: i64,
    coeffs: Vec<(i64, Q)>
}

pub fn print_gap(z: &Number) -> String {
    let mut str_list: Vec<String> = vec![];
    for (exp, coeff) in &z.coeffs {
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
    pub fn new(order: i64, coeffs: &Vec<(i64, Q)>) -> Number {
        Number {
            order: order,
            coeffs: coeffs.clone(),
        }
    }

    pub fn increase_order_to(z: &mut Self, new_order: i64) {
        for (&mut exp, _) in &mut z.coeffs {
            exp = new_order * exp / z.order;
        }
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

fn add_single(coeffs: &mut Vec<(i64, Q)>, exp: i64, coeff: &Q, sign: Sign) {
    // TODO: this is horrifically bad for efficiency, please replace
    let mut exp_index = -1;
    let sign = if sign == Sign::Plus { 1 } else { -1 };

    for i in 0..coeffs.len() {
        // Found the exponent, no need to shift stuff around, we need to
        // actually add.
        if coeffs[i].0 == exp {
            exp_index = i;
            let exp_coeff = &mut coeffs[i].1;
            exp_coeff.add_assign(sign * coeff);
            return;
        }

        // exp didn't appear in coeffs.
        // This is the first bigger exponent we've seen, so we know
        // exp belongs at index i.
        if coeffs[i].0 > exp {
            exp_index = i;
            coeffs.insert(i, (exp, sign * coeff.clone()));
            return;
        }
    }

    // exp is larger than any exponent in coeffs right now
    if exp_index == -1 {
        coeffs.push((exp, sign*coeff));
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

impl Arbitrary for Number {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        random_cyclotomic(g, 100, 101)
    }
}

field_axiom_tests!(Number);
