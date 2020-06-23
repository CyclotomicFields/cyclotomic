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
use std::convert::TryInto;
use std::thread::current;

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
    num::integer::gcd(*x, *y) == 1
}

fn phi(n: &u64) -> u64 {
    let mut count = 0;
    for k in 1..*n {
        if are_coprime(n, &k) {
            println!("n = {:?} coprime to k = {:?}", n, k);
            count += 1;
        }
    }
    count
}

fn get_same_coeff(z: &Number) -> Option<Q> {
    let coeffs = z.coeffs.clone().into_iter().map(|(exp, coeff)| { coeff });
    let nonzero_coeffs: HashSet<Q> = coeffs.filter(|q| { *q != Q::zero() }).collect();

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
            None => current_gcd = Some(*exp),
            Some(gcd) => current_gcd = Some(num::integer::gcd(gcd, *exp))
        }
    }

    if current_gcd.is_none() {
        // just zero appeared, so this is a rational, we can reduce all the way
        // to E(2)
        if saw_zero {
            let coeff = z.coeffs.get(&0).unwrap();
            return Number::new(&2, &vec![(0, coeff.clone())].into_iter().collect());
        }

        // not even zero appeared, so this is just zero
        return zero_order(&2);
    }

    // otherwise gcd is Some and nonzero, so we can divide through

    println!("gcd is {:?}", current_gcd.unwrap());

    let mut new_order = z.order.clone() / current_gcd.unwrap();
    let mut new_coeffs = z.coeffs.clone().into_iter().map(|p| { (p.0 / current_gcd.unwrap(), p.1) }).collect();
    let res = Number::new(&new_order, &new_coeffs);

    println!("reduced {:?} to {:?}", z, res);
    res
}

fn try_rational(z: &Number) -> Option<Q> {
    println!("orig z = {:?}", z);
    let base_z = convert_to_base(&try_reduce(&convert_to_base(z)));
    println!("base_z = {:?}", base_z);

    let n = base_z.order as i64;
    let phi_n = phi(&(n as u64));
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

    let is_squarefree = n_div_powers.clone().into_iter().all(|(factor, power)| { power < 2 });
    println!("is_squarefree = {:?}", is_squarefree);

    let num_primes = n_div_powers.clone().into_iter().filter(|(factor, power)| { power > &(0 as i64) }).count();
    println!("num_primes = {:?}", num_primes);

    let num_nonzero_terms = base_z.coeffs.clone().into_iter().filter(|(exp, coeff)| { *coeff != Q::zero() }).count();
    println!("num_nonzero_terms = {:?}", num_nonzero_terms);

    let same_coeff = get_same_coeff(&base_z);
    println!("same_coeff = {:?}", same_coeff);

    match same_coeff {
        Some(coeff) => if num_nonzero_terms == phi_n as usize && is_squarefree {
            println!("it's rational!");
            Some(coeff.mul(Z::from(i64::pow(-1, num_primes.try_into().unwrap()))))
        } else {
            None
        }
        None => None
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

    let n = z.order.clone() as i64;

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
            let c = result.coeffs.get(&(i as u64)).unwrap_or(&zero).clone();
            if c == Q::zero() {
                continue;
            }

            let orig = Number::e(n as u64, i as u64).scalar_mul(&c);

            println!("i = {:?}", i);

            if set
                .clone()
                .into_iter()
                .any(|x| math_mod(&x, &q) == math_mod(&i, &q))
            {
                println!("i = {:?} not ok, rewriting", i);

                // delete from result, we'll add it back later rewritten in the
                // basis
                result.coeffs.remove(&(i as u64));

                // use the relation to rewrite this root, these are the
                // exponents on the right hand side
                let exps: Vec<i64> = (1..*p).map(|k| math_mod(&(k * n / p + i), &n)).collect();

                let mut rhs = zero_order(&(n as u64));
                for k in exps {
                    rhs.coeffs.insert(k as u64, -c.clone());
                }
                println!("rewriting {:?} = {:?}", orig, rhs);
                result = result + rhs;
            } else {
                // just because i is ok with regard to p, it might be bad
                // with regard to a different, later prime. We leave it alone.
                println!("i = {:?} ok w.r.t p = {:?}, leaving alone", i, p);
            }
        }
    }

    result
}

/// Expresses a cyclotomic as an element of the smallest cyclotomic field
/// possible.
fn reduce_to_min_subfield(z: &Number) -> Number {
    // TODO: write this. but it's only needed for efficiency, not correctness
    z.clone()
}

fn apply_automorphism(z: &Number, i: &u64) -> Number {
    let mut result = zero_order(&z.order);

    for (exp, coeff) in z.coeffs.clone() {
        result.coeffs.insert(math_mod(&((exp * (*i)) as i64), &(z.order as i64)) as u64, coeff);
    }

    result
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
        let n = &z.order;

        // The Galois group $G = \text{Aut}(L/K)$ has order $\phi(n)$. The
        // elements are the automorphisms $\zeta_n \mapsto \zeta_n^i$ for all
        // $1 \leq i \leq n-1$ coprime to $n$.

        // This is the product except for the term for $t = \id_L$.
        let mut x = one_order(&z.order);

        for i in 2..*n {
            if are_coprime(&i, n) {
                x = x.clone() * apply_automorphism(z, &i);
            }
        }
        println!("x = {:?}", x);

        // The full product:
        let q_cyc = z.clone() * x.clone();
        println!("q_cyc = {:?}", q_cyc);

        // q_cyc is rational, so let's extract the rational bit (all of it).
        // We need to do some tricks to make it rational (it might be in a
        // bit of a weird form).

        // How do you tell if a cyclotomic is really rational?
        let q: Q = try_rational(&q_cyc).unwrap();
        println!("q = {:?}", q);

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
    fn e(n: u64, k: u64) -> Self {
        Number::new(
            &n,
            &[(k, Q::from_integer(Z::one()))].iter().cloned().collect(),
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
    fn try_rational_test() {
        let q = Q::new(Z::from(10), Z::from(3));
        println!("q = {:?}", q);
        let z = one_order(&10).scalar_mul(&q);
        println!("z = {:?}", z);
        let z_rat = try_rational(&z);
        println!("z_rat = {:?}", z_rat);
    }

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

    #[test]
    fn mul_has_inverses_specific() {
        let z = Number::e(10, 3);
        let prod = convert_to_base(&(z.inv() * z.clone()));
        println!("prod = {:?}", prod);
        println!("is_one = {:?}", prod == one_order(&z.order));
    }

    quickcheck! {
    fn mul_has_inverses(z: Number) -> bool {
        // TODO: skip if z is zero, no inverse
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

    #[test]
    fn mul_distribute_specific() {
        let e = Number::e;

        let x = e(3, 1) + e(3, 2);
        let y = e(4, 2);
        let z = e(3, 1) + e(3, 2);

        let lhs = x.clone() * (y.clone() + z.clone());
        let rhs = x.clone() * y.clone() + x.clone() * z.clone();

        let real_lhs = convert_to_base(&lhs);
        let real_rhs = convert_to_base(&rhs);

        println!("orig lhs = {:?}", lhs);
        println!("basis lhs = {:?}", real_lhs);

        println!("orig rhs = {:?}", rhs);
        println!("basis rhs = {:?}", real_rhs);
    }

    #[test]
    fn convert_to_base_specific() {
        let e = Number::e;
        let x = e(15, 3) + e(15, 5);
        println!("E(15)^3 + E(15)^5 = {:?}", convert_to_base(&x));
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
            let orders: Vec<u64> = vec![3, 5, 6, 7, 10];
            let order = orders[g.gen_range(0, orders.len())];
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
