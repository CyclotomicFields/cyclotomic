// Utilities for computing with Galois groups of field extensions over Q
// Probably only Abelian ones, i.e. cyclotomic fields.

use crate::fields::big_sparse::Number;
use crate::fields::util::*;
use crate::fields::CyclotomicFieldElement;
use crate::fields::big_sparse::Exponent;
use std::convert::TryInto;

pub fn apply_automorphism(z: &Number, i: &Exponent) -> Number {
    // TODO: really support bigints
    let mut result = Number::zero_order((&z.order).try_into().unwrap());

    for (exp, coeff) in z.coeffs.clone() {
        result.coeffs.insert(math_mod_big(&(&exp * i).into(), &z.order), coeff);
    }

    result
}
