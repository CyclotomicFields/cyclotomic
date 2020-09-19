// Utilities for computing with Galois groups of field extensions over Q
// Probably only Abelian ones, i.e. cyclotomic fields.

use crate::fields::exponent::Exponent;
use crate::fields::rational::Rational;
use crate::fields::sparse::Number;
use crate::fields::util::*;
use crate::fields::CyclotomicFieldElement;
use std::convert::TryInto;

pub fn apply_automorphism<E: Exponent, Q: Rational>(z: &Number<E, Q>, i: &E) -> Number<E, Q> {
    let mut result = Number::zero_order(&z.order);

    for (exp, coeff) in z.coeffs.clone() {
        result.coeffs.insert(
            Exponent::math_mod(&(exp.clone() * i.clone()).into(), &z.order),
            coeff,
        );
    }

    result
}
