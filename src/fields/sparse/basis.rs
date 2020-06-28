// Functions that require specific knowledge of which basis is used for the
// field. In this case, we use the Zumbroich basis, described in GAP
// documentation and probably some papers.

use super::num::Zero;
use crate::fields::sparse::*;
use crate::fields::*;
use std::collections::HashSet;
use std::convert::TryInto;
use std::ops::Mul;
use fnv::FnvHashMap;

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

    let new_order = z.order.clone() / current_gcd.unwrap() as i64;

    let mut new_coeffs = FnvHashMap::default();
    for (exp, coeff) in &z.coeffs {
        if coeff.is_zero() {
            continue;
        }
        let new_exp = *exp / current_gcd.unwrap() as i64;
        new_coeffs.insert(new_exp, coeff.clone());
    }
    let res = Number::new(new_order, &new_coeffs);

    res
}

pub fn try_rational(z: &Number) -> Option<Q> {
    let base_z = convert_to_base(&try_reduce(&convert_to_base(z)));

    // if order is 2, that's already a rational!
    if base_z.order == 2 {
        return Some(base_z.coeffs.get(&0).unwrap_or(&Q::zero()).clone());
    }

    let n = base_z.order;
    let phi_n = phi(n);

    let mut n_divisors: Vec<i64> = divisors::get_divisors(n as u64)
        .into_iter()
        .map(|x| x as i64)
        .collect();

    if n_divisors.len() == 0 {
        n_divisors.push(n);
    }

    let n_div_powers = &count_powers(&n, &n_divisors);

    let is_squarefree = n_div_powers
        .into_iter()
        .all(|(_, power)| *power < 2);

    let num_primes = n_div_powers
        .into_iter()
        .filter(|(factor, power)| *power > 0)
        .count();

    let num_nonzero_terms = base_z
        .coeffs
        .clone()
        .into_iter()
        .filter(|(exp, coeff)| *coeff != Q::zero())
        .count();

    let same_coeff = get_same_coeff(&base_z);

    match same_coeff {
        Some(coeff) => {
            if num_nonzero_terms == phi_n as usize && is_squarefree {
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

    let n = z.order;

    // currently z, will by the end still be equal to z but will be written in
    // the Zumbroich basis
    let mut result = z.clone();

    for (exp, coeff) in &z.coeffs {
        if coeff == &Q::zero() {
            result.coeffs.remove(exp);
        }
    }

    let n_divisors: Vec<i64> = divisors::get_divisors(n as u64)
        .into_iter()
        .map(|x| x as i64)
        .collect();
    let mut n_div_powers = count_powers(&n, &n_divisors);

    // if it has no divisors smaller than itself, it's prime
    if n_div_powers.is_empty() {
        n_div_powers.push((n, 1));
    }

    // TODO: rewrite this to be more functional and readable?
    for (p, power) in &n_div_powers {
        // the maximal power of p that divides n
        let q: i64 = p.pow(*power as u32);

        // i is in this set (mod q) iff it is not a basis element
        let start_bad = if *p == 2 { q / 2 } else { -(q / p - 1) / 2 };
        let end_bad = if *p == 2 { q - 1 } else { (q / p - 1) / 2 };
        let mut bad_exponents = Vec::<i64>::new();
        bad_exponents.reserve((end_bad - start_bad + 1) as usize);
        for x in start_bad..=end_bad {
            bad_exponents.push(n / q * x);
        }

        for i in 0..n {
            // if there isn't even a term for i, why consider it?
            let coeff = {
                let maybe_coeff = result.coeffs.get(&i);

                match maybe_coeff {
                    None => continue,
                    Some(rational) => {
                        if rational.is_zero() {
                            continue;
                        }
                    }
                }

                maybe_coeff.unwrap().clone()
            };

            for bad_exp in &bad_exponents {
                if math_mod(&bad_exp, &q) == math_mod(&i, &q) {
                    // delete from result, we'll add it back later, rewritten
                    result.coeffs.remove(&i);

                    // use the relation to rewrite this root
                    for k in 1..*p {
                        let new_exp = math_mod(&(k * n / p + i), &n);
                        let maybe_existing_coeff = result.coeffs.get(&new_exp).clone();

                        match maybe_existing_coeff {
                            None => {
                                result.coeffs.insert(new_exp, -coeff.clone());
                            }
                            Some(existing_coeff) => {
                                if existing_coeff.is_zero() {
                                    result.coeffs.insert(new_exp, -coeff.clone());
                                } else {
                                    result
                                        .coeffs
                                        .insert(new_exp, existing_coeff - coeff.clone());
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    result
}
