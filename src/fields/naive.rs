extern crate num;

use self::num::{One, Zero};
use crate::fields::{CyclotomicFieldElement, Exponent, FieldElement, Q, Z};
use num::traits::Inv;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use std::cmp::Eq;
use std::collections::HashMap;
use std::collections::HashSet;
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
    pub fn new(order: &Exponent, coeffs: &HashMap<Exponent, Q>) -> Number {
        Number {
            order: order.clone(),
            coeffs: coeffs.clone(),
        }
    }
    pub fn increase_order_to(z: &mut Self, new_order: &u64) -> () {
        let mut new_coeffs = HashMap::new();
        for (exp, coeff) in z.coeffs.clone() {
            new_coeffs.insert(new_order * exp / z.order, coeff);
        }
        z.order = new_order.clone();
        z.coeffs = new_coeffs;
    }

    pub fn match_orders(z1: &mut Number, z2: &mut Number) -> () {
        let new_order = num::integer::lcm(z1.order, z2.order);
        Number::increase_order_to(z1, &new_order);
        Number::increase_order_to(z2, &new_order);
        assert_eq!(z1.order, z2.order);
    }
}

/// Represents zero iff all coeffs are zero.
impl Zero for Number {
    fn zero() -> Self {
        Number::new(&3, &HashMap::new())
    }

    fn set_zero(&mut self) {
        self.coeffs = HashMap::new();
    }

    fn is_zero(&self) -> bool {
        self.coeffs.values().all(|coeff| coeff.is_zero())
    }
}

// Creates a zero already written as an element of $\mathbb{Q}(\zeta_n)$.
fn zero_order(n: &u64) -> Number {
    Number::new(n, &HashMap::new())
}

impl One for Number {
    fn one() -> Self {
        // seems like a good default - the smallest nontrivial cyclotomic field
        one_order(&3)
    }
}

// Creates a one already written as an element of $\mathbb{Q}(\zeta_n)$.
fn one_order(n: &u64) -> Number {
    let mut coeffs = HashMap::new();
    for i in 1..n.clone() {
        coeffs.insert(i, -Q::one());
    }
    Number::new(n, &coeffs)
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
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(&mut z1, &mut z2);

        // We will never need to reduce here, you can't add low powers of
        // $\zeta_n$ and get higher powers. Higher powers do not exist.
        let mut coeffs: HashMap<u64, Q> = HashMap::new();
        for (exp, coeff) in z1.coeffs.into_iter().chain(z2.coeffs) {
            match coeffs.clone().get(&exp) {
                Some(existing_coeff) => coeffs.insert(exp, coeff + existing_coeff),
                None => coeffs.insert(exp, coeff),
            };
        }

        let result = Number::new(&z1.order, &coeffs);
        result
    }
}

impl Mul for Number {
    type Output = Number;

    /// Multiplies term by term, not bothering to do anything interesting.
    fn mul(self, rhs: Self) -> Self::Output {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(&mut z1, &mut z2);

        let mut result = zero_order(&z1.order);

        // This order is almost certainly not optimal. But you know, whatever.
        result.order = z1.order;
        for (exp1, coeff1) in z1.coeffs.clone() {
            for (exp2, coeff2) in z2.coeffs.clone() {
                let new_exp = (exp1 + exp2) % z1.order.clone();
                let new_coeff = coeff1.clone() * coeff2.clone();

                // Special case: if the new exponent would be 0, since 1 is not
                // a basis element, we have to use the fact that:
                // $1 = -\sum_{i=1}^{p-1} \zeta_n^i$ to rewrite the new constant
                // term in our basis.
                if new_exp != 0 {
                    match result.coeffs.clone().get(&new_exp) {
                        Some(existing_coeff) => {
                            result.coeffs.insert(new_exp, new_coeff + existing_coeff)
                        }
                        None => result.coeffs.insert(new_exp, new_coeff),
                    };
                } else {
                    for i in 1..result.order.clone() {
                        match result.coeffs.clone().get(&i) {
                            Some(existing_coeff) => {
                                result.coeffs.insert(i, existing_coeff - new_coeff.clone())
                            }
                            None => result.coeffs.insert(i, -new_coeff.clone()),
                        };
                    }
                }
            }
        }
        result
    }
}

fn are_coprime(x: &u64, y: &u64) -> bool {
    let x_divs = divisors::get_divisors(*x);
    let y_divs = divisors::get_divisors(*y);

    for div in x_divs {
        if y_divs.contains(&div) {
            return false;
        }
    }

    return true;
}

fn phi(n: &u64) -> u64 {
    let mut count = 0;
    for k in 1..*n {
        if are_coprime(n, &k) {
            count += 1;
        }
    }
    count
}

// TODO: smaller functions!!! more tests!!!
fn try_rational(z: &Number) -> Option<Q> {
    let p = &z.order;

    // if $z.order = p$ is prime, then there is only 1 nontrivial orbit, it's
    // $\zeta_p^i$ for $1 \leq i \leq p-1$, it's of size $p-1$. For z to be
    // rational, it must be $p + q(\sum_{i=1}^{p-1} \zeta_p^i)$ where $p, q$ are
    // rational.

    // The coefficients of the $\zeta_p^i$ with $i > 0$
    let mut root_coeffs = z.coeffs.clone();
    root_coeffs.remove(&0);

    println!("root_coeffs = {:?}", root_coeffs);

    // all nonzero powers must appear and must have the same coeff
    let mut all_same_coeff = root_coeffs.get(&1).unwrap_or(&Q::zero()).clone();

    println!("first coeff = {:?}", all_same_coeff.to_string());

    for i in 2..*p {
        let coeff = root_coeffs.get(&i);
        println!("coeff = {:?}", coeff);
        match coeff {
            // This root doesn't appear, and we saw another root with nonzero
            // coefficient, bad! z can't be rational!
            None => {
                if !all_same_coeff.clone().is_zero() {
                    return None;
                }
            }

            // This root does appear (with possibly zero coeff), so we need it
            // to match the other coeffs we've seen
            Some(q) => {
                if all_same_coeff != *q {
                    return None;
                }
            }
        }

        // if we get to here, the coeff matched all others we've seen
        all_same_coeff = coeff.unwrap_or(&Q::zero()).clone();
    }

    println!("all coeffs matched = {:?}", all_same_coeff);

    // if we got here it means all the coeffs matched
    // the -1 comes from the fact that adding all of the nontrivial roots
    // gives -1
    return Some(all_same_coeff * -Q::one());
}

fn counts(v: &Vec<u64>) -> HashMap<i64, u64> {
    let mut v_counts: HashMap<i64, u64> = HashMap::new();

    for elem in v {
        let current_count = v_counts.get(&(*elem as i64)).unwrap_or(&0);
        v_counts.insert(*elem as i64, current_count + 1);
    }

    v_counts
}

fn math_mod(x: &i64, n: &i64) -> i64 {
    (x % n + n) % n
}

fn is_in_basis(i: &i64, n: &i64, n_div_powers: &HashMap<i64, u64>) -> bool {
    for (p, power) in n_div_powers {
        // the maximal power of p that divides n
        let q: i64 = p ^ (&(*power as i64));

        // i is in this set (mod q) iff it is not a basis element
        let set: HashSet<i64> = if *p == 2 {
            (q / 2..q - 1).map(|x| *n / q * x).collect()
        } else {
            (-(q / p - 1) / 2..(q / p - 1) / 2)
                .map(|x| *n / q * x)
                .collect()
        };

        if set.into_iter().any(|x| math_mod(&x, &q) == math_mod(i, &q)) {
            return false;
        }
    }
    true
}

// THIS IMPLEMENTATION OF THE BASIS STUFF IS PURE SHIT
// TODO: improve it by reading gap's cyclotom.c?

/// Returns a vector of the i such that the $\zeta_n^i$ form the GAP basis
/// for $K$.
fn gap_basis(n: &u64) -> Vec<u64> {
    let n_divisors = divisors::get_divisors(*n);
    let n_div_powers = counts(&n_divisors);

    let mut basis_elems = vec![];

    // TODO: get rid of nested loops and casts
    for i in 0..=*n {
        if is_in_basis(&(i as i64), &(*n as i64), &n_div_powers) {
            basis_elems.push(i);
        }
    }

    basis_elems
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
        let z = self;
        println!("z = {:?}", z);

        // Let $L = \mathbb{Q}(\zeta_n), K = \mathbb{Q}$.
        // Then $L/K$ is a degree $\phi(n)$ extension.

        // I can't be bothered to implement $\phi$ so let's say the order is
        // prime so $\phi(n) = n-1$ lmao
        // TODO: don't be lazy!

        let n = &z.order;

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
        let g = |z: &Number| {
            // $g$ doubles all powers of $\zeta_n$.
            // Note that it's a field automorphism so permutes the roots,
            // we don't have to worry about collisions or anything like that.
            let mut result = zero_order(&z.order);

            for (exp, coeff) in &z.coeffs {
                result.coeffs.insert((exp * 2) % n, coeff.clone());
            }

            result
        };

        // $\prod_{t \in G} t(z)$ is rational, and we have $G$ cyclic, of
        // order $p-1$ so:
        // $\prod{i=0}^{p-2} g^i(z) = q \in K$.

        // Applied f times times to z
        fn apply_times(f: impl Fn(&Number) -> Number, times: &u64, z: &Number) -> Number {
            let mut result = z.clone();
            for i in 0..*times {
                result = f(&result);
            }
            result
        }

        // This is the product except for the term for $t = \id_L$.
        let mut x = one_order(&z.order);
        for i in 1..n - 1 {
            let term = apply_times(g, &i, &z);
            x = x * term;
        }

        // The full product:
        let q_cyc = z.clone() * x.clone();
        println!("q_cyc = {:?}", q_cyc);

        // q_cyc is rational, so let's extract the rational bit (all of it).
        // We need to do some tricks to make it rational (it might be in a
        // bit of a weird form).
        let q: Q = try_rational(&q_cyc).unwrap();
        println!("q = {:?}", q);

        // Now the fun part:
        // $\prod{i=0}^{p-2} g^i(z) = q$ so $z (q^{-1} \prod{i=1}^{p-2} g^i(z)) = 1$.
        // So $z^{-1} = q^{-1} \prod{i=1}^{p-1} g^i(z) = x/q$.
        let z_inv = x.scalar_mul(&q.inv());

        println!("z_inv = {:?}", z_inv);

        z_inv
    }

    /// Also could probably be better.
    fn inv_mut(&mut self) -> Self {
        self.inv()
    }
}

impl CyclotomicFieldElement for Number {
    fn e(n: &u64, k: &u64) -> Self {
        Number::new(
            n,
            &[(k.clone(), Q::from_integer(Z::one()))]
                .iter()
                .cloned()
                .collect(),
        )
    }

    fn scalar_mul(&self, scalar: &Q) -> Self {
        let mut result = self.clone();
        for (_, coeff) in result.coeffs.iter_mut() {
            *coeff *= scalar.clone();
        }
        result
    }
}

// These are precisely the field axioms.
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
        z.clone() + zero_order(&z.order) == z.clone()
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
        let same = z.clone() * one_order(&z.order);
        println!("same = {:?}", same);
        same == z.clone()
    }
    }

    quickcheck! {
    fn add_has_inverses(z: Number) -> bool {
        let minus_one = Q::new(Z::from(-1), Z::from(1));
        z.clone() + z.clone().scalar_mul(&minus_one) == zero_order(&z.order)
    }
    }

    quickcheck! {
    fn zero_kills_all(z: Number) -> bool {
        zero_order(&z.order) * z.clone() == zero_order(&z.order)
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
        let prod = z.inv() * z.clone();
        println!("prod = {:?}", prod);
        prod == one_order(&z.order)
    }
    }

    quickcheck! {
    fn mul_distributes_over_add(x: Number, y: Number, z: Number) -> bool {
        x.clone() * (y.clone() + z.clone()) == x.clone() * y.clone() + x.clone() * z.clone()
    }
    }

    // These are some non field axiom, implementation-specific tests
    #[test]
    fn test_gap_basis() {
        // This is the example from GAP
        let basis_45: HashSet<u64> = gap_basis(&45).into_iter().collect();
        let expected: HashSet<u64> = vec![
            1, 2, 3, 6, 7, 8, 11, 12, 16, 17, 19, 21, 24, 26, 28, 29, 33, 34, 37, 38, 39, 42, 43,
            44,
        ]
        .into_iter()
        .collect();
    }

    // Just a hack to get arbitrary coefficients
    #[derive(Clone)]
    struct QArb(Q);

    impl Arbitrary for QArb {
        fn arbitrary<G>(g: &mut G) -> Self
        where
            G: Gen,
        {
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
            let order = 5; //orders[g.gen_range(0, orders.len())];
            let num_terms: u64 = g.gen_range(1, 5);
            let mut result = zero_order(&order);

            for _ in 1..=num_terms {
                let exp: u64 = g.gen_range(1, order);
                let QArb(coeff) = QArb::arbitrary(g);
                result.coeffs.insert(exp, coeff);
            }

            result
        }
    }
}
