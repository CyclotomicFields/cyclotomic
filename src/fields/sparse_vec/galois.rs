// Utilities for computing with Galois groups of field extensions over Q
// Probably only Abelian ones, i.e. cyclotomic fields.

use crate::fields::sparse_vec::{math_mod, Number};
use crate::fields::CyclotomicFieldElement;

pub fn apply_automorphism(z: &Number, i: i64) -> Number {
    let mut result = Number::zero_order(z.order);

    for (exp, coeff) in z.coeffs.clone() {
        result.coeffs.push((math_mod(&(exp * i), &z.order), coeff));
    }

    // this is O(nlogn), surely not optimal! TODO: improve!
    result.coeffs.sort_by(|(exp1, _), (exp2, _)| { exp1.cmp(exp2) });

    result
}
