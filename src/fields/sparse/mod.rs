extern crate num;

use self::num::{One, Zero};
use num::traits::Inv;
use crate::fields::{CyclotomicFieldElement, FieldElement, Q, Z};
use crate::fields::{AdditiveGroup, MultiplicativeGroup};
use quickcheck::{Arbitrary, Gen};
use rand::Rng;
use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::TryInto;
use std::fmt;
use std::vec::Vec;
use std::ops::Mul;

pub mod add;
pub mod mul;

/// Represents a polynomial in the `order`th root of unity.
///
/// Simplest possible choice of data structure - a hash map. We assume only that
/// each term has an exponent less than the order. So the only real term is the
/// one for exponent zero, for example.
#[derive(Clone)]
pub struct Number {
    order: i64,
    coeffs: HashMap<i64, Q>,
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
    pub fn new(order: i64, coeffs: &HashMap<i64, Q>) -> Number {
        Number {
            order: order,
            coeffs: coeffs.clone(),
        }
    }

    pub fn increase_order_to(z: &mut Self, new_order: i64) {
        let mut new_coeffs = HashMap::new();
        for (exp, coeff) in &z.coeffs {
            new_coeffs.insert(new_order * exp.clone() / z.order.clone(), coeff.clone());
        }
        z.order = new_order.clone();
        z.coeffs = new_coeffs;
    }

    pub fn match_orders(z1: &mut Number, z2: &mut Number) {
        let new_order = num::integer::lcm(z1.order, z2.order);
        Number::increase_order_to(z1, new_order);
        Number::increase_order_to(z2, new_order);
        assert_eq!(z1.order, z2.order);
    }
}

// Creates a zero already written as an element of $\mathbb{Q}(\zeta_n)$.
fn zero_order(n: i64) -> Number {
    Number::new(n, &HashMap::new())
}

// Creates a one already written as an element of $\mathbb{Q}(\zeta_n)$.
fn one_order(n: i64) -> Number {
    let mut coeffs = HashMap::new();
    for i in 1..n.clone() {
        coeffs.insert(i, -Q::one());
    }
    Number::new(n, &coeffs)
}

fn are_coprime(x: i64, y: i64) -> bool {
    num::integer::gcd(x as u64, y as u64) == 1
}

fn phi(n: i64) -> i64 {
    let mut count = 0;
    for k in 1..n {
        if are_coprime(n, k) {
            println!("n = {:?} coprime to k = {:?}", n, k);
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

// Tries to reduce to a possibly smaller cyclotomic field
fn try_reduce(z: &Number) -> Number {
    let mut current_gcd: Option<u64> = None;
    let mut saw_zero = false;

    for (exp, coeff) in &z.coeffs {
        // this term doesn't really appear
        if *coeff == Q::zero() {
            continue;
        }

        // 0 will mess up the gcd calculation
        if *exp == 0 {
            saw_zero = true;
            continue;
        }

        match current_gcd {
            None => current_gcd = Some(*exp as u64),
            Some(gcd) => current_gcd = Some(num::integer::gcd(gcd, *exp as u64)),
        }
    }

    if current_gcd.is_none() {
        // just zero appeared, so this is a rational, we can reduce all the way
        // to E(2)
        if saw_zero {
            let coeff = z.coeffs.get(&0).unwrap();
            return Number::new(2, &vec![(0, coeff.clone())].into_iter().collect());
        }

        // not even zero appeared, so this is just zero
        return zero_order(2);
    }

    // otherwise gcd is Some and nonzero, so we can divide through

    println!("gcd is {:?}", current_gcd.unwrap());

    let new_order = z.order.clone() / current_gcd.unwrap() as i64;

    let mut new_coeffs = HashMap::new();
    for (exp, coeff) in &z.coeffs {
        if coeff.is_zero() {
            continue;
        }
        let new_exp = *exp / current_gcd.unwrap() as i64;
        println!("new_exp = {:?}, coeff = {:?}", new_exp, coeff.clone());
        new_coeffs.insert(new_exp, coeff.clone());
    }
    let res = Number::new(new_order, &new_coeffs);

    println!("new_coeffs = {:?}", new_coeffs);
    println!("reduced {:?} to {:?}", z, res);
    res
}

fn try_rational(z: &Number) -> Option<Q> {
    println!("orig z = {:?}", z);
    let base_z = convert_to_base(&try_reduce(&convert_to_base(z)));
    println!("base_z = {:?}", base_z);

    // if order is 2, that's already a rational!
    if base_z.order == 2 {
        return Some(base_z.coeffs.get(&0).unwrap_or(&Q::zero()).clone());
    }

    let n = base_z.order;
    let phi_n = phi(n);
    println!("phi(n) = {:?}", phi_n);

    let mut n_divisors: Vec<i64> = divisors::get_divisors(n as u64)
        .into_iter()
        .map(|x| x as i64)
        .collect();

    if n_divisors.len() == 0 {
        n_divisors.push(n);
    }

    let n_div_powers = count_powers(&n, &n_divisors);
    println!("n divisors = {:?}", n_div_powers);

    let is_squarefree = n_div_powers
        .clone()
        .into_iter()
        .all(|(factor, power)| power < 2);
    println!("is_squarefree = {:?}", is_squarefree);

    let num_primes = n_div_powers
        .clone()
        .into_iter()
        .filter(|(factor, power)| *power > 0)
        .count();
    println!("num_primes = {:?}", num_primes);

    let num_nonzero_terms = base_z
        .coeffs
        .clone()
        .into_iter()
        .filter(|(exp, coeff)| *coeff != Q::zero())
        .count();
    println!("num_nonzero_terms = {:?}", num_nonzero_terms);

    let same_coeff = get_same_coeff(&base_z);
    println!("same_coeff = {:?}", same_coeff);

    match same_coeff {
        Some(coeff) => {
            if num_nonzero_terms == phi_n as usize && is_squarefree {
                println!("it's rational!");
                Some(coeff.mul(Z::from(i64::pow(-1, num_primes.try_into().unwrap()))))
            } else {
                None
            }
        }
        None => None,
    }
}

fn math_mod(x: &i64, n: &i64) -> i64 {
    let res = (x % n + n) % n;
    println!("{:?} mod {:?} is {:?}", x, n, res);
    res
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

/// Writes a cyclotomic in the GAP basis (keeps field the same). This is needed
/// for elements to have a unique representation, e.g. to check equality
/// The rewriting rule is taken from GAP but the actual code isn't.
fn convert_to_base(z: &Number) -> Number {
    // Currently z is expressed as a sum of some $e_n^i$.
    // We need to eliminate the $i$ such that either:
    // $i \in n/q [-(q/p-1)/2 .. (q/p-1)/2] mod q
    // for some odd prime p dividing n with q = p^k with k maximal s.t. q | n
    // or $i \in n/q [q/2 .. q-1]$ mod q where q = 2^k with k maximal s.t. q | n

    // Suppose $i$ satisfies one of the conditions because of a prime p.
    // To get rid of $c e_n^i$ (c rational), we write the identity
    // $1+e_p+e_p^2+..+e_p^{p-1}=0$ in $n$th roots, so it's:
    // $0=1+e_n^{n/p}+e_n^{2n/p}+..+e_n^{(p-1)n/p}$
    // and so: $0=e_n^i+e_n^{n/p+i}+e_n^{2n/p+i}+..+e_n^{(p-1)n/p+i}$
    // then finally: $c e_n^i = -c(e_n^{n/p+i}+e_n^{2n/p+i}+..+e_n^{(p-1)n/p+i})$
    // All of the powers on the right are in the basis, too, but it's
    // not trivial, there's a bit of an argument there, maybe a proof.

    // Note: this is written for readability and is too functional for its
    // own good. The machine code generated might not be very good.

    let n = z.order.clone();

    // currently z, will by the end still be equal to z but will be written in
    // the Zumbroich basis
    let mut result = z.clone();

    let n_divisors: Vec<i64> = divisors::get_divisors(n as u64)
        .into_iter()
        .map(|x| x as i64)
        .collect();
    let mut n_div_powers = count_powers(&n, &n_divisors);
    println!("divisors = {:?}", n_divisors);

    // if it has no divisors smaller than itself, it's prime
    if n_div_powers.is_empty() {
        n_div_powers.push((n, 1));
    }

    println!("n is {:?}", n);
    println!("n divisors: {:?}", n_div_powers);

    // TODO: rewrite this to be more functional and readable?
    for (p, power) in &n_div_powers {
        // the maximal power of p that divides n
        let q: i64 = p.pow(*power as u32);
        println!("q = {:?}^{:?} = {:?}", p, power, q);

        // i is in this set (mod q) iff it is not a basis element
        let set: HashSet<i64> = if *p == 2 {
            (q / 2..q - 1 + 1).map(|x| n / q * x).collect()
        } else {
            (-(q / p - 1) / 2..(q / p - 1) / 2 + 1)
                .map(|x| n / q * x)
                .collect()
        };

        println!("bad i are {:?}", set);

        let zero = Q::zero();

        for i in 0..n.clone() {
            // if there isn't even a term for i, why consider it?
            let c = result.coeffs.get(&i).unwrap_or(&zero).clone();
            if c == Q::zero() {
                continue;
            }

            let orig = Number::e(n, i).scalar_mul(&c).clone();

            println!("i = {:?}", i);

            if set
                .clone()
                .into_iter()
                .any(|x| math_mod(&x, &q) == math_mod(&i, &q))
            {
                println!("i = {:?} not ok, rewriting", i);

                // delete from result, we'll add it back later rewritten in the
                // basis
                result.coeffs.remove(&i);

                // use the relation to rewrite this root, these are the
                // exponents on the right hand side
                let exps: Vec<i64> = (1..*p).map(|k| math_mod(&(k * n / p + i), &n)).collect();

                let mut rhs = zero_order(n.clone());
                for k in exps {
                    rhs.coeffs.insert(k, -c.clone());
                }
                println!("rewriting {:?} = {:?}", orig, rhs);
                result.add(&mut rhs);
            } else {
                // just because i is ok with regard to p, it might be bad
                // with regard to a different, later prime. We leave it alone.
                println!("i = {:?} ok w.r.t p = {:?}, leaving alone", i, p);
            }
        }
    }

    result
}

fn apply_automorphism(z: &Number, i: i64) -> Number {
    let mut result = zero_order(z.order);

    for (exp, coeff) in z.coeffs.clone() {
        result.coeffs.insert(
            math_mod(&(exp * i), &z.order),
            coeff,
        );
    }

    result
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
}

// These are precisely the field axioms.
// TODO: write these tests for the trait, then instantiate somehow? can rust do that?
#[cfg(test)]
mod tests {
    use super::*;

    quickcheck! {
    fn zero_is_add_identity(z: Number) -> bool {
        z.clone().add(&mut zero_order(z.order.clone())).eq(&mut z.clone())
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
        let mut same = z.clone().mul(&mut one_order(z.order.clone())).clone();
        println!("same = {:?}", same);
        same.eq(&mut z.clone())
    }
    }

    quickcheck! {
    fn add_has_inverses(z: Number) -> bool {
        z.clone().add(z.clone().add_invert()).eq(&mut zero_order(z.order))
    }
    }

    quickcheck! {
    fn zero_kills_all(z: Number) -> bool {
        zero_order(z.order.clone()).mul(&mut z.clone()).eq(&mut zero_order(z.order))
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
    fn mul_has_inverses(z: Number) -> bool {
        // TODO: skip if z is zero, no inverse
        let mut prod = z.clone().mul_invert().mul(&mut z.clone()).clone();
        println!("prod = {:?}", prod);
        prod.eq(&mut one_order(z.order))
    }
    }

    quickcheck! {
    fn mul_distributes_over_add(x: Number, y: Number, z: Number) -> bool {
        x.clone().mul(y.clone().add(&mut z.clone())).eq(x.clone().mul(&mut y.clone()).add(x.clone().mul(&mut z.clone())))
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
            let p: i64 = g.gen_range(1, 10);
            let q: i64 = g.gen_range(1, 10);
            QArb(Q::new(Z::from(p), Z::from(q)))
        }
    }

    impl Arbitrary for Number {
        fn arbitrary<G>(g: &mut G) -> Self
        where
            G: Gen,
        {
            let orders: Vec<i64> = vec![3, 4, 5, 6, 7, 8, 9, 10];
            let order = orders[g.gen_range(0, orders.len())];
            let num_terms: u64 = g.gen_range(1, 5);
            let mut result = zero_order(order.clone());

            for _ in 1..=num_terms {
                let exp: i64 = g.gen_range(1, order);
                let QArb(coeff) = QArb::arbitrary(g);
                result.coeffs.insert(exp, coeff);
            }

            result
        }
    }
}
