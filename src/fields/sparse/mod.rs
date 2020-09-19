extern crate num;
extern crate rustc_hash;

use self::num::{One, Zero};

use crate::fields::rational::Rational;
use crate::fields::util::*;
use crate::fields::MultiplicativeGroupElement;
use crate::fields::{CyclotomicFieldElement, FieldElement, Z};
use basis::convert_to_base;
use num::traits::Inv;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::{AddAssign, Mul, SubAssign};
use std::vec::Vec;

#[macro_use]
use crate::fields::*;
use self::rustc_hash::FxHashMap;
use crate::fields::exponent::Exponent;

pub mod add;
pub mod basis;
pub mod galois;
pub mod mul;

// TODO: how to make this FxHashMap for i64 and HashMap for Z?
type ExpCoeffMap<E, Q> = FxHashMap<E, Q>;

/// Represents a polynomial in the `order`th root of unity.
#[derive(Clone)]
pub struct Number<E: Exponent = i64, Q: Rational = rug::Rational> {
    order: E,
    pub coeffs: ExpCoeffMap<E, Q>,
}

pub fn print_gap<E: Exponent, Q: Rational>(z: &Number<E, Q>) -> String {
    let mut str_list: Vec<String> = vec![];
    let mut exp = E::from(0);
    while &exp != &z.order {
        let zero = Q::zero();
        let coeff = z.coeffs.get(&exp).unwrap_or(&zero);
        if !coeff.is_zero() {
            str_list.push(String::from(
                format!("{} * E({})^{}", coeff, z.order, exp).as_str(),
            ))
        }
        exp = exp + E::from(1);
    }
    "(".to_string() + &str_list.join(" + ") + ")"
}

impl<E, Q> fmt::Debug for Number<E, Q>
where
    E: Exponent,
    Q: Rational,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Number ({})", print_gap(self))
    }
}

impl<E, Q> Number<E, Q>
where
    E: Exponent,
    Q: Rational,
{
    pub fn new(order: &E, coeffs: &ExpCoeffMap<E, Q>) -> Number<E, Q> {
        Number {
            order: order.clone(),
            coeffs: coeffs.clone(),
        }
    }

    pub fn increase_order_to(z: &mut Self, new_order: &E) {
        let mut new_coeffs = ExpCoeffMap::default();
        for (exp, coeff) in &z.coeffs {
            new_coeffs.insert(
                new_order.clone() * exp.clone() / z.order.clone(),
                coeff.clone(),
            );
        }
        z.order = new_order.clone();
        z.coeffs = new_coeffs;
    }

    pub fn match_orders(z1: &mut Number<E, Q>, z2: &mut Number<E, Q>) {
        if z1.order == z2.order {
            return;
        }
        let new_order: E = Exponent::lcm(&z1.order, &z2.order);
        Number::<E, Q>::increase_order_to(z1, &new_order);
        Number::<E, Q>::increase_order_to(z2, &new_order);
    }
}

fn get_same_coeff<E: Exponent, Q: Rational>(z: &Number<E, Q>) -> Option<Q> {
    let coeffs = z.coeffs.clone().into_iter().map(|(_exp, coeff)| coeff);
    let nonzero_coeffs: HashSet<Q> = coeffs.filter(|q| !q.is_zero()).collect();

    if nonzero_coeffs.len() == 0 {
        // all coeffs are zero
        Some(Q::zero())
    } else if nonzero_coeffs.len() == 1 {
        Some(nonzero_coeffs.iter().last().unwrap().clone())
    } else {
        None
    }
}

fn add_single<E: Exponent, Q: Rational>(
    coeffs: &mut ExpCoeffMap<E, Q>,
    exp: &E,
    coeff: &Q,
    sign: Sign,
) {
    let maybe_existing_coeff = coeffs.get_mut(exp);
    match maybe_existing_coeff {
        None => {
            if sign == Sign::Plus {
                coeffs.insert(exp.clone(), coeff.clone());
            } else {
                let mut neg = coeff.clone();
                neg.add_invert();
                coeffs.insert(exp.clone(), neg);
            }
        }
        Some(existing_coeff) => {
            // TODO: find a way to get rid of coeff.clone() here, it's not needed
            if sign == Sign::Plus {
                existing_coeff.add(&mut coeff.clone());
            } else {
                existing_coeff.add(coeff.clone().add_invert());
            }
        }
    }
}

pub fn is_zero<E: Exponent, Q: Rational>(z: &Number<E, Q>) -> bool {
    for (_, coeff) in &z.coeffs {
        if !coeff.is_zero() {
            return false;
        }
    }
    true
}

impl<E, Q> FieldElement for Number<E, Q>
where
    E: Exponent,
    Q: Rational,
{
    fn eq(&mut self, other: &mut Self) -> bool {
        let mut za = self.clone();
        let mut zb = other.clone();
        Number::<E, Q>::match_orders(&mut za, &mut zb);
        let z1 = convert_to_base(&za);
        let z2 = convert_to_base(&zb);

        // Now that we've matched the orders, z1 and z2 are expressed as
        // elements in the same field so are the same iff each nonzero term is
        // the same.
        fn has_diff<E: Exponent, Q: Rational>(left: &Number<E, Q>, right: &Number<E, Q>) -> bool {
            for (exp_left, coeff_left) in &left.coeffs {
                match right.coeffs.get(&exp_left) {
                    None => {
                        if !coeff_left.is_zero() {
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

impl<E, Q> CyclotomicFieldElement<E, Q> for Number<E, Q>
where
    E: Exponent,
    Q: Rational,
{
    fn e(n: &E, k: &E) -> Self {
        Number::<E, Q>::new(n, &[(k.clone(), Q::from((1, 1)))].iter().cloned().collect())
    }

    fn scalar_mul(&mut self, scalar: &Q) -> &mut Self {
        let mut result = self.clone();
        for (_, coeff) in result.coeffs.iter_mut() {
            // TODO: REMOVE THIS CLONE
            coeff.mul(&mut scalar.clone());
        }
        *self = result;
        self
    }

    fn zero_order(n: &E) -> Number<E, Q> {
        Number::<E, Q>::new(&n, &ExpCoeffMap::<E, Q>::default())
    }

    fn one_order(n: &E) -> Number<E, Q> {
        let mut coeffs = ExpCoeffMap::default();
        let mut i = E::from(1);
        while i != *n {
            coeffs.insert(i.clone(), Q::from((-1, 1)));
            i = i + E::from(1);
        }
        Number::<E, Q>::new(n, &coeffs)
    }

    fn complex_conjugate(&self) -> Self {
        let mut new_coeffs = ExpCoeffMap::default();

        for (exp, coeff) in &self.coeffs {
            if *exp == E::from(0) {
                new_coeffs.insert(E::from(0), coeff.clone());
            } else {
                new_coeffs.insert(self.order.clone() - exp.clone(), coeff.clone());
            }
        }

        Self::new(&self.order, &new_coeffs)
    }
}

pub fn random_rational<G, Q: Rational>(g: &mut G) -> Q
where
    G: rand::RngCore,
{
    let p: i64 = g.gen_range(1, 10);
    let q: u64 = g.gen_range(1, 10);
    Q::from((p, q))
}

pub fn random_cyclotomic<G, E, Q>(g: &mut G, min_order: i64, max_order: i64) -> Number<E, Q>
where
    G: rand::RngCore,
    E: Exponent,
    Q: Rational,
{
    let order = g.gen_range(min_order, max_order);
    let num_terms: u64 = g.gen_range(1, 5);
    let mut result = Number::<E, Q>::zero_order(&E::from(order.clone()));

    for _ in 1..=num_terms {
        let exp: i64 = g.gen_range(1, order);
        let coeff = random_rational(g);
        result.coeffs.insert(E::from(exp), coeff);
    }

    result
}

impl<E: 'static, Q: 'static> Arbitrary for Number<E, Q>
where
    E: Exponent,
    Q: Rational,
{
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        random_cyclotomic::<G, E, Q>(g, 2, 50)
    }
}

type Number_i64 = Number<i64, rug::Rational>;
type Number_Z = Number<Z, rug::Rational>;

#[cfg(test)]
mod i64_tests {
    use super::*;
    field_axiom_tests!(Number_i64);
}

#[cfg(test)]
mod Z_tests {
    use super::*;
    field_axiom_tests!(Number_Z);
}
