// Utilities for computing with Galois groups of field extensions over Q
// Probably only Abelian ones, i.e. cyclotomic fields.

use crate::fields::sparse::{math_mod, Number};
use crate::fields::CyclotomicFieldElement;

pub fn apply_automorphism(z: &Number, i: i64) -> Number {
    let mut result = Number::zero_order(z.order);

    for (exp, coeff) in z.coeffs.clone() {
        result.coeffs.insert(math_mod(&(exp * i), &z.order), coeff);
    }

    result
}
