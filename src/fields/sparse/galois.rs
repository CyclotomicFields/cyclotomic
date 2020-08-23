// Utilities for computing with Galois groups of field extensions over Q
// Probably only Abelian ones, i.e. cyclotomic fields.

use crate::fields::sparse::Number;
use crate::fields::util::*;
use crate::fields::CyclotomicFieldElement;
use std::convert::TryInto;
use crate::fields::exponent::Exponent;

pub fn apply_automorphism<E: Exponent>(z: &Number<E>, i: &E) -> Number<E> {
    let mut result = Number::zero_order(z.order.clone());

    for (exp, coeff) in z.coeffs.clone() {
        result.coeffs.insert(Exponent::math_mod(&(&exp * i).into(), &z.order), coeff);
    }

    result
}
