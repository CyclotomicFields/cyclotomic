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
use rustc_hash::FxHashMap;
use std::collections::HashSet;
use std::fmt;
use std::ops::{AddAssign, Mul, SubAssign};
use std::vec::Vec;

#[macro_use]
use crate::fields::*;

pub mod add;
pub mod basis;
pub mod galois;
pub mod mul;

pub type Exponent = Z;
pub type ExpCoeffMap = FxHashMap<Exponent, Q>;

/// Represents a polynomial in the `order`th root of unity.
#[derive(Clone)]
pub struct Number {
    order: Exponent,
    coeffs: ExpCoeffMap,
}

pub fn print_gap(z: &Number) -> String {
    let mut str_list: Vec<String> = vec![];
    let mut exp = Exponent::from(0);
    while &exp != &z.order {
        let zero = Q::from(0).clone();
        let coeff = z.coeffs.get(&exp).unwrap_or(&zero);
        if *coeff != 0 {
            str_list.push(String::from(
                format!("{} * E({})^{}", coeff, z.order, exp).as_str(),
            ))
        }
        exp += 1;
    }
    "(".to_string() + &str_list.join(" + ") + ")"
}

impl fmt::Debug for Number {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Number ({})", print_gap(self))
    }
}

impl Number {
    pub fn new(order: &Exponent, coeffs: &ExpCoeffMap) -> Number {
        Number {
            order: order.clone(),
            coeffs: coeffs.clone(),
        }
    }

    pub fn increase_order_to(z: &mut Self, new_order: &Exponent) {
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
        let new_order: Exponent = z1.order.lcm_ref(&z2.order).into();
        Number::increase_order_to(z1, &new_order);
        Number::increase_order_to(z2, &new_order);
    }
}

fn get_same_coeff(z: &Number) -> Option<Q> {
    let coeffs = z.coeffs.clone().into_iter().map(|(_exp, coeff)| coeff);
    let nonzero_coeffs: HashSet<Q> = coeffs.filter(|q| *q != 0).collect();

    if nonzero_coeffs.len() == 0 {
        // all coeffs are zero
        Some(Q::from(0))
    } else if nonzero_coeffs.len() == 1 {
        Some(nonzero_coeffs.iter().last().unwrap().clone())
    } else {
        None
    }
}

fn add_single(coeffs: &mut ExpCoeffMap, exp: &Exponent, coeff: &Q, sign: Sign) {
    let maybe_existing_coeff = coeffs.get_mut(exp);
    match maybe_existing_coeff {
        None => {
            if sign == Sign::Plus {
                coeffs.insert(exp.clone(), coeff.clone());
            } else {
                coeffs.insert(exp.clone(), -coeff.clone());
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
        if *coeff != 0 {
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

        // Now that we've matched the orders, z1 and z2 are expressed as
        // elements in the same field so are the same iff each nonzero term is
        // the same.
        fn has_diff(left: &Number, right: &Number) -> bool {
            for (exp_left, coeff_left) in &left.coeffs {
                match right.coeffs.get(&exp_left) {
                    None => {
                        if *coeff_left != 0 {
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
            &Exponent::from(n),
            &[(k.into(), Q::from(1))].iter().cloned().collect(),
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
        Number::new(&n.into(), &ExpCoeffMap::default())
    }

    fn one_order(n: i64) -> Number {
        let mut coeffs = ExpCoeffMap::default();
        for i in 1..n {
            coeffs.insert(i.into(), Q::from(-1));
        }
        Number::new(&n.into(), &coeffs)
    }
}

pub fn random_rational<G>(g: &mut G) -> Q
where
    G: rand::RngCore,
{
    let p: i64 = g.gen_range(1, 10);
    let q: i64 = g.gen_range(1, 10);
    Q::from((p, q))
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
        result.coeffs.insert(Exponent::from(exp), coeff);
    }

    result
}

impl Arbitrary for Number {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        random_cyclotomic(g, 2, 50)
    }
}

//field_axiom_tests!(Number);
