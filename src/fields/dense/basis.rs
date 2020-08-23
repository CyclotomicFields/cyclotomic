// Functions that require specific knowledge of which basis is used for the
// field. In this case, we use the Zumbroich basis, described in GAP
// documentation and probably some papers.

use super::num::Zero;
use crate::fields::dense::*;

use std::convert::TryInto;
use std::ops::Mul;
use crate::fields::util::Sign;
use crate::fields::exponent::Exponent;
use crate::fields::CyclotomicFieldElement;

// Tries to reduce to a possibly smaller cyclotomic field
pub fn try_reduce(z: &mut Number) {
    let mut current_gcd: Option<u64> = None;
    let mut saw_exp_zero = false;
    let mut num_nonzero_terms = 0;
    let mut coeffs_are_equal = true;
    let mut last_nonzero_coeff: Option<Q> = None;

    for exp in 0..z.order {
        let coeff = z.coeffs[exp as usize].clone();
        // this term doesn't really appear
        if coeff == 0 {
            continue;
        }

        num_nonzero_terms += 1;

        if coeffs_are_equal {
            match &last_nonzero_coeff {
                None => last_nonzero_coeff = Some(coeff.clone()),
                Some(seen_coeff) => {
                    if *seen_coeff != coeff {
                        coeffs_are_equal = false;
                    }
                }
            }
        }

        // 0 will mess up the gcd calculation
        if exp == 0 {
            saw_exp_zero = true;
            continue;
        }

        match current_gcd {
            None => current_gcd = Some(exp as u64),
            Some(gcd) => current_gcd = Some(num::integer::gcd(gcd, exp as u64)),
        }
    }

    if current_gcd.is_none() {
        // if the current gcd was never set, then either 0 is the only exponent
        // or there are no exponents - rational in both cases.
        z.order = 1;

        if saw_exp_zero {
            let coeff = z.coeffs[0].clone();
            *z = Number::zero_order(&z.order);
            z.coeffs[0] = coeff;
        } else {
            *z = Number::zero_order(&z.order);
        }

        return;
    }

    // otherwise gcd is Some and nonzero, so we can divide through
    let gcd = current_gcd.unwrap() as i64;

    // no point dividing through by 1
    if gcd != 1 {
        let new_order = z.order / gcd;
        z.order = new_order;

        // don't need to reduce exp=0 since it would just reduce to exp=0
        for exp in 1..new_order {
            z.coeffs[exp as usize] = z.coeffs[(gcd * exp) as usize].clone();
        }

        z.coeffs.truncate(new_order as usize);
    }

    // now all exponents are coprime with the new order
    let n = z.order;
    let phi_n = Exponent::phi(&n);

    let n_div_powers = Exponent::factorise(&n);
    let is_squarefree = n_div_powers.clone().into_iter().all(|(_, power)| power < 2);
    let num_primes = n_div_powers
        .into_iter()
        .filter(|(_factor, power)| *power > 0)
        .count();

    // if this is the case, it's rational
    if num_nonzero_terms == phi_n && coeffs_are_equal && is_squarefree {
        z.order = 1;
        z.coeffs.clear();
        let new_coeff = last_nonzero_coeff
            .unwrap()
            .mul(Z::from(i64::pow(-1, num_primes.try_into().unwrap())));
        z.coeffs.push(new_coeff);
        return;
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
    let mut n_div_powers = Exponent::factorise(&n);


    for (p, power) in &n_div_powers {
        // the maximal power of p that divides n
        let q: i64 = p.pow(*power as u32);

        // i is in this set (mod q) iff it is not a basis element
        let start_bad = if *p == 2 { q / 2 } else { -(q / p - 1) / 2 };
        let end_bad = if *p == 2 { q - 1 } else { (q / p - 1) / 2 };

        for bad_exp_raw in start_bad..=end_bad {
            let bad_exp = Exponent::math_mod(&(n / q * bad_exp_raw), &q);
            // We want to remove the i that are equal to bad_exp mod q.
            // These are exactly the i such that i = bad_exp + aq for some a.
            // We also need only check 0 <= i <= n-1, which means we only need
            // to check -bad_exp/q - 1 <= a <= (n-1-bad_exp)/q + 1.
            // The -1 and +1 are so that even if the division isn't perfect,
            // then we still check the full range of a we need to check.
            for a in -bad_exp / q - 1..=((n - 1 - bad_exp) / q + 1) {
                let i: i64 = bad_exp + a * q;
                // if there isn't even a term for i, no need to convert it
                let coeff = result.coeffs[(Exponent::math_mod(&i, &n)) as usize].clone();

                if coeff == 0 {
                    continue;
                }

                // if we got here, i has a nonzero term so must be rewritten
                result.coeffs[Exponent::math_mod(&i, &n) as usize] = Q::from(0);
                for k in 1..*p {
                    let new_exp = Exponent::math_mod(&(k * n / p + i), &n);
                    add_single(&mut result.coeffs, new_exp, &coeff, Sign::Minus);
                }
            }
        }
    }

    result
}
