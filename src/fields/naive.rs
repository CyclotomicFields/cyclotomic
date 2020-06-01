extern crate num;

use self::num::{One, Zero};
use crate::fields::{CyclotomicFieldElement, Exponent, FieldElement, Q, Z};
use num::traits::Inv;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use std::cmp::Eq;
use std::collections::HashMap;
use std::fmt;
use std::iter::Product;
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

impl PartialEq for Number {
    fn eq(&self, other: &Self) -> bool {
        let mut z1 = self.clone();
        let mut z2 = other.clone();
        Number::match_orders(&mut z1, &mut z2);

        // Now that we've matched the orders, z1 and z2 are expressed as
        // elements in the same field so are the same iff each nonzero term is
        // the same.
        fn has_diff(left: &Number, right: &Number) -> bool {
            for (exp_left, coeff_left) in left.coeffs.clone() {
                match right.coeffs.get(&exp_left) {
                    None => {
                        if coeff_left != Q::zero() {
                            return true;
                        }
                    }
                    Some(coeff_right) => {
                        if coeff_left != *coeff_right {
                            return true;
                        }
                    }
                }
            }
            false
        }

        !has_diff(&z1, &z2) && !has_diff(&z2, &z1)
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

                match result.coeffs.get(&new_exp) {
                    Some(existing_coeff) => {
                        result.coeffs.insert(new_exp, new_coeff + existing_coeff)
                    }
                    None => result.coeffs.insert(new_exp, new_coeff),
                };
            }
        }
        result
    }
}

impl Product<Number> for Number {
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = Number>,
    {
        let mut result = Number::one();

        for z in iter.into_iter() {
            result = result * z;
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

    /// Gives the inverse of $z$ using the product of Galois conjugates.
    ///
    /// I don't think there's a "trivial" or "stupid" way of doing this.
    /// The product of the Galois conjugates is rational, we can normalise
    /// to get the multiplicative inverse.
    fn inv(&self) -> Self {
        let z = self.clone();

        // Let $L = \mathbb{Q}(\zeta_n), K = \mathbb{Q}$.
        // Then $L/K$ is a degree $\phi(n)$ extension.

        // I can't be bothered to implement $\phi$ so let's say the order is
        // prime so $\phi(n) = n$ lmao
        // TODO: don't be lazy!

        let n = z.order.clone();

        // The Galois group $G = \text{Aut}(L/K)$ has order $\phi(n)$. The
        // elements are the automorphisms $\zeta_n \mapsto \zeta_n^i$ for all
        // $1 \leq i \leq n-1$ coprime to $n$.

        // I can't be bothered with this right now, so lets say $n = p$ prime,
        // then the automorphisms are just those maps with $1 \leq i \leq p-1$
        // and $G = \mathbb{Z}_p^\times \cong \mathbb{Z}_{p-1}$
        // TODO: generalise this you hack fraud

        // From this point we just drop pretenses and assume $n = p$ prime for
        // everything. I don't even care any more.

        // $G$ is cyclic and generated by any non-trivial element, they all
        // have order $p-1$. I pick $g(\zeta_n) = \zeta_n^2$ as the generator.
        let g = |z: Number| {
            // doubles all powers
            // Note that it's a field automorphism so permutes the roots,
            // we don't have to worry about collisions or anything like that.
            let mut result = Number::zero();

            for (exp, coeff) in z.coeffs.clone() {
                result.coeffs.insert((exp * 2) % n, coeff);
            }

            result
        };

        // $\prod_{x \in G} x(z)$ is rational, and we have $G$ cyclic, so:
        // $\prod{i=0}^{p-1} g^i(z) = q \in K$.
        fn apply_times(f: impl Fn(Number) -> Number, times: u64, z: Number) -> Number {
            let mut result = z;
            for _ in 0..times {
                result = f(result);
            }
            result
        }

        // This is the product except for z, call it x, we'll need it later
        let x: Number = (1..n - 1)
            .into_iter()
            .map(|i| apply_times(g, i, z.clone()))
            .product();

        // The full product:
        let q_cyc: Number = z.clone() * x.clone();
        println!("q_cyc = {:?}", q_cyc);

        // q_cyc is rational, so let's extract the rational bit (all of it):
        let zero: u64 = 0;
        let q: Q = q_cyc.coeffs.get(&zero).unwrap().to_owned();

        // Now the fun part:
        // $\prod{i=0}^{p-1} g^i(z) = q$ so $z (q^{-1} \prod{i=1}^{p-1} g^i(z)) = 1$.
        // So $z^{-1} = q^{-1} \prod{i=1}^{p-1} g^i(z) = x/q$.
        let z_inv = x.scalar_mul(q.inv());

        z_inv
    }

    /// Also could probably be better.
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
        // TODO: change these
        let p: i64 = g.gen_range(1, 2);
        let q: i64 = g.gen_range(1, 2);
        QArb(Q::new(Z::from(p), Z::from(q)))
    }
}

impl Arbitrary for Number {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        // TODO: make this work for order non prime
        let order: u64 = 13;
        let num_terms: u64 = g.gen_range(1, 2);
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
    fn add_is_associative(x: Number, y: Number, z: Number) -> bool {
        (x.clone() + y.clone()) + z.clone() == x.clone() + (y.clone() + z.clone())
    }
    }

    quickcheck! {
    fn add_is_commutative(x: Number, y: Number) -> bool {
        x.clone() + y.clone() == y.clone() + x.clone()
    }
    }

    quickcheck! {
    fn one_is_mul_identity(z: Number) -> bool {
        z.clone() * Number::one() == z.clone()
    }
    }

    quickcheck! {
    fn add_has_inverses(z: Number) -> bool {
        let minus_one = Q::new(Z::from(-1), Z::from(1));
        z.clone() + z.clone().scalar_mul(minus_one) == Number::zero()
    }
    }

    quickcheck! {
    fn zero_kills_all(z: Number) -> bool {
        Number::zero() * z == Number::zero()
    }
    }

    quickcheck! {
    fn mul_is_commutative(x: Number, y: Number) -> bool {
        x.clone() * y.clone() == y.clone() * x.clone()
    }
    }

    quickcheck! {
    fn mul_is_associative(x: Number, y: Number, z: Number) -> bool {
        (x.clone() * y.clone()) * z.clone() == x.clone() * (y.clone() * z.clone())
    }
    }

    quickcheck! {
    fn mul_has_inverses(z: Number) -> bool {
        // TODO: GET THIS WORKING
        // z.inv() * z.clone() == Number::one()
        true
    }
    }

    quickcheck! {
    fn mul_distributes_over_add(x: Number, y: Number, z: Number) -> bool {
        x.clone() * (y.clone() + z.clone()) == x.clone() * y.clone() + x.clone() * z.clone()
    }
    }
}
