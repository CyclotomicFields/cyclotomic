// Functions that require specific knowledge of which basis is used for the
// field. In this case, we use the Zumbroich basis, described in GAP
// documentation and probably some papers.

use crate::fields::*;
use crate::fields::sparse::*;
use super::num::Zero;
use std::collections::{HashMap, HashSet};
use std::ops::Mul;
use std::convert::TryInto;

// Tries to reduce to a possibly smaller cyclotomic field
pub fn try_reduce(z: &Number) -> Number {
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
        return Number::zero_order(2);
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

pub fn try_rational(z: &Number) -> Option<Q> {
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

/// Writes a cyclotomic in the GAP basis (keeps field the same). This is needed
/// for elements to have a unique representation, e.g. to check equality
/// The rewriting rule is taken from GAP but the actual code isn't.
pub fn convert_to_base(z: &Number) -> Number {
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

                let mut rhs = Number::zero_order(n.clone());
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
