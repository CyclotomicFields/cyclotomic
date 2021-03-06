// Utilities for computing with Galois groups of field extensions over Q
// Probably only Abelian ones, i.e. cyclotomic fields.

use crate::fields::dense::Number;
use crate::fields::util::*;
use crate::fields::CyclotomicFieldElement;
use crate::fields::exponent::Exponent;

pub fn apply_automorphism(z: &Number, i: i64) -> Number {
    let mut result = Number::zero_order(&z.order);

    for exp in 0..z.order {
        result.coeffs[Exponent::math_mod(&(exp * i), &z.order) as usize] = z.coeffs[exp as usize].clone();
    }

    result
}
