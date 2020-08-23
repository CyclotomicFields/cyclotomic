// Functions that require specific knowledge of which basis is used for the
// field. In this case, we use the Zumbroich basis, described in GAP
// documentation and probably some papers.

use super::num::Zero;
use crate::fields::sparse::*;

use rug::ops::Pow;
use crate::fields::util::Sign;
use std::convert::TryInto;
use std::ops::Mul;
use crate::fields::exponent::Exponent;

// Tries to reduce to a possibly smaller cyclotomic field
pub fn try_reduce<E: Exponent>(z: &mut Number<E>) {
    let mut current_gcd: Option<E> = None;
    let mut saw_exp_zero = false;
    let mut coeffs_are_equal = true;
    let mut last_nonzero_coeff: Option<Q> = None;
    let mut num_nonzero_terms = E::from(0);

    for (exp, coeff) in &z.coeffs {
        // this term doesn't really appear
        if *coeff == 0 {
            continue;
        }

        num_nonzero_terms += 1;

        if coeffs_are_equal {
            match &last_nonzero_coeff {
                None => last_nonzero_coeff = Some(coeff.clone()),
                Some(seen_coeff) => {
                    if seen_coeff != coeff {
                        coeffs_are_equal = false;
                    }
                }
            }
        }

        // 0 will mess up the gcd calculation
        if *exp == 0 {
            saw_exp_zero = true;
            continue;
        }

        match current_gcd {
            None => current_gcd = Some(exp.clone()),
            Some(gcd) => current_gcd = Some(gcd.gcd_ref(exp).into()),
        }
    }

    if current_gcd.is_none() {
        // if the current gcd was never set, then either 0 is the only exponent
        // or there are no exponents - rational in both cases.
        z.order = E::from(1);

        if saw_exp_zero {
            let coeff = z.coeffs.get(&E::from(0)).unwrap().clone();
            z.coeffs.clear();
            z.coeffs.insert(E::from(0), coeff);
        } else {
            z.coeffs.clear();
        }

        return;
    }

    // otherwise gcd is Some and nonzero, so we can divide through
    let gcd = current_gcd.unwrap();

    // no point dividing through by 1
    if gcd != 1 {
        let new_order: E = (&z.order / &gcd).into();
        z.order = new_order.clone();

        // don't need to reduce exp=0 since it would just reduce to exp=0
        let mut exp = E::from(1);
        while &exp != &new_order {
            match z.coeffs.get(&(&gcd * &exp).into()) {
                None => (), // there is no term to rewrite
                Some(coeff) => {
                    z.coeffs.insert(exp.clone(), coeff.clone());
                    z.coeffs.remove(&(&gcd * &exp).into());
                }
            }
            exp += 1;
        }
    }

    // now all exponents are coprime with the new order
    let n = &z.order;
    let phi_n = Exponent::phi(n);

    let n_div_powers = &Exponent::factorise(n);
    let is_squarefree = n_div_powers.into_iter().all(|(_, power)| *power < 2);
    let num_primes = n_div_powers
        .into_iter()
        .filter(|(_factor, power)| *power > 0)
        .count();

    // if this is the case, it's rational
    if num_nonzero_terms == phi_n && coeffs_are_equal && is_squarefree {
        z.order = E::from(1);
        z.coeffs.clear();
        let new_coeff = last_nonzero_coeff
            .unwrap()
            .mul(Z::from(i64::pow(-1, num_primes.try_into().unwrap())));
        z.coeffs.insert(E::from(0), new_coeff);
        return;
    }
}

/// Writes a cyclotomic in the GAP basis (keeps field the same). This is needed
/// for elements to have a unique representation, e.g. to check equality
/// The rewriting rule is taken from GAP but the actual code isn't.
pub fn convert_to_base<E>(z: &Number<E>) -> Number<E> where E: Exponent {
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

    let n = &z.order;

    // currently z, will by the end still be equal to z but will be written in
    // the Zumbroich basis
    let mut result = z.clone();

    for (exp, coeff) in &z.coeffs {
        if *coeff == 0 {
            result.coeffs.remove(exp);
        }
    }

    let mut n_div_powers = Exponent::factorise(n);

    for (p, power) in &n_div_powers {
        // the maximal power of p that divides n
        let q: E = p.pow(*power).into();

        // i is in this set (mod q) iff it is not a basis element
        let start_bad: E = if *p == 2 {
            &q / E::from(2)
        } else {
            let div: E = (&q / p).into();
            let minus1: E = (div - E::from(1)).into();
            let neg: E = (-minus1).into();
            (neg / E::from(2)).into()
        };
        let end_bad: E = if *p == 2 {
            &q - E::from(1)
        } else {
            ((&q / p).into(): E - 1) / 2
        };

        let mut bad_exp_raw = start_bad;
        while &bad_exp_raw <= &end_bad {
            let bad_exp = Exponent::math_mod(&((n / &q).into(): E * &bad_exp_raw), &q);
            // We want to remove the i that are equal to bad_exp mod q.
            // These are exactly the i such that i = bad_exp + aq for some a.
            // We also need only check 0 <= i <= n-1, which means we only need
            // to check -bad_exp/q - 1 <= a <= (n-1-bad_exp)/q + 1.
            // The -1 and +1 are so that even if the division isn't perfect,
            // then we still check the full range of a we need to check.
            let start_check = (-&bad_exp).into(): E / &q - E::from(1);
            let end_check =
                (((n - E::from(1)).into(): E - &bad_exp) / &q + E::from(1));
            let mut a = start_check;
            while a <= end_check {
                let i = (&bad_exp + &a * &q).into();
                // if there isn't even a term for i, no need to convert it
                let coeff = {
                    let maybe_coeff = result.coeffs.get(&i);

                    match maybe_coeff {
                        None => {
                            a += 1;
                            continue;
                        }
                        Some(rational) => {
                            if *rational == 0 {
                                a += 1;
                                continue;
                            }
                        }
                    }
                    maybe_coeff.unwrap().clone()
                };

                // if we got here, i has a nonzero term so must be rewritten
                result.coeffs.remove(&i);
                let mut k = Z::from(1);
                while &k != p {
                    let new_exp = Exponent::math_mod(&((&k * n).into(): E / p + &i), &n);
                    add_single(&mut result.coeffs, &new_exp, &coeff, Sign::Minus);
                    k += 1;
                }
                a += 1;
            }
            bad_exp_raw += 1;
        }
    }

    result
}
