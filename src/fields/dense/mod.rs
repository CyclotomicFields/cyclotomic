extern crate num;
extern crate rustc_hash;

use self::num::{One, Zero};
use crate::fields::util::*;
use crate::fields::MultiplicativeGroupElement;
use crate::fields::{CyclotomicFieldElement, FieldElement, Q, Z};
use basis::convert_to_base;
use num::traits::Inv;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use std::fmt;
use std::ops::{AddAssign, Mul, SubAssign};
use std::vec::Vec;

#[macro_use]
use crate::fields::*;

#[macro_use]
use std::collections::HashSet;

pub mod add;
pub mod basis;
pub mod galois;
pub mod mul;

/// Represents a polynomial in the `order`th root of unity.
#[derive(Clone)]
pub struct Number {
    order: i64,
    pub coeffs: Vec<Q>,
}

pub fn print_gap(z: &Number) -> String {
    let mut str_list: Vec<String> = vec![];
    for exp in 0..z.order {
        let coeff = &z.coeffs[exp as usize];
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
    pub fn new(order: i64, coeffs: &Vec<Q>) -> Number {
        Number {
            order: order,
            coeffs: coeffs.clone(),
        }
    }

    pub fn increase_order_to(z: &mut Self, new_order: i64) {
        let mut new_coeffs = Vec::with_capacity(new_order as usize);
        // not the most cache-friendly. TODO: improve?
        for _new_exp in 0..new_order {
            new_coeffs.push(Q::zero());
        }
        for old_exp in 0..z.order {
            new_coeffs[(new_order * old_exp / z.order) as usize] =
                z.coeffs[old_exp as usize].clone();
        }
        z.order = new_order;
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

fn get_same_coeff(z: &Number) -> Option<Q> {
    let nonzero_coeffs: HashSet<Q> = z
        .coeffs
        .clone()
        .into_iter()
        .filter(|q| !q.is_zero())
        .collect();

    if nonzero_coeffs.len() == 0 {
        // all coeffs are zero
        Some(Q::zero())
    } else if nonzero_coeffs.len() == 1 {
        Some(nonzero_coeffs.iter().last().unwrap().clone())
    } else {
        None
    }
}

fn add_single(coeffs: &mut Vec<Q>, exp: i64, coeff: &Q, sign: Sign) {
    let existing_coeff = &mut coeffs[exp as usize];
    if sign == Sign::Plus {
        *existing_coeff += coeff;
    } else {
        *existing_coeff -= coeff;
    }
}

pub fn is_zero(z: &Number) -> bool {
    for coeff in &z.coeffs {
        if !coeff.is_zero() {
            return false;
        }
    }
    true
}

impl FieldElement for Number {
    fn eq(&mut self, other: &mut Self) -> bool {
        let mut za = self.clone();
        let mut zb = other.clone();
        Number::match_orders(&mut za, &mut zb);
        let z1 = convert_to_base(&za);
        let z2 = convert_to_base(&zb);

        for i in 0..z1.order {
            if z1.coeffs[i as usize] != z2.coeffs[i as usize] {
                return false;
            }
        }
        true
    }
}

impl CyclotomicFieldElement for Number {
    fn e(n: i64, k: i64) -> Self {
        let mut coeffs = vec![Q::zero(); n as usize];
        coeffs[k as usize] = Q::from_integer(Z::one());
        Number::new(n, &coeffs)
    }

    fn scalar_mul(&mut self, scalar: &Q) -> &mut Self {
        for coeff in &mut self.coeffs {
            *coeff *= scalar;
        }
        self
    }

    fn zero_order(n: i64) -> Number {
        Number::new(n, &vec![Q::zero(); n as usize])
    }

    fn one_order(n: i64) -> Number {
        let mut coeffs = vec![Q::zero(); n as usize];
        coeffs[0] = Q::one();
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
        result.coeffs[exp as usize] = coeff;
    }

    result
}

impl Arbitrary for Number {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        random_cyclotomic(g, 2, 10)
    }
}

field_axiom_tests!(Number);
