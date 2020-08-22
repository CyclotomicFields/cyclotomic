pub type Z = rug::Integer;
pub type Q = rug::Rational;

use num::Integer;
use std::collections::HashMap;

pub trait AdditiveGroupElement {
    /// Adds z to self in place, so self = self + z
    fn add(&mut self, z: &mut Self) -> &mut Self;

    /// Negates z in in place
    fn add_invert(&mut self) -> &mut Self;
}

pub trait MultiplicativeGroupElement {
    /// Multiplies self by z in place, so self = self * z
    fn mul(&mut self, z: &mut Self) -> &mut Self;

    /// Inverts self in place.
    fn mul_invert(&mut self) -> &mut Self;
}

/// Provides operations for fields. Expected to satisfy the field axioms.
pub trait FieldElement: AdditiveGroupElement + MultiplicativeGroupElement {
    /// Equality, but can shuffle coefficients around and simplify expressions
    /// for greater efficiency.
    fn eq(&mut self, other: &mut Self) -> bool;
}

/// Provides convenience functions specific to cyclotomic fields.
pub trait CyclotomicFieldElement<Exponent>: FieldElement + Clone
    where Exponent: From<i64>
{
    /// Returns $\zeta_n$^k.
    fn e(n: &Exponent, k: &Exponent) -> Self;

    /// Multiplies in-place by scalar. Recall that $\mathbb{Q}(\zeta_n)$ is a
    /// $\mathbb{Q}$-vector space.
    fn scalar_mul(&mut self, scalar: &Q) -> &mut Self;

    /// Gives zero expressed as an element of $\mathbb{Q}(\zeta_n)$
    fn zero_order(n: &Exponent) -> Self;

    /// Gives one expressed as an element of $\mathbb{Q}(\zeta_n)$
    fn one_order(n: &Exponent) -> Self;
}

/// Possible data structure for a CyclotomicFieldElement, useful as a common
/// interchange format for different concrete CyclotomicFieldElement types.
#[derive(Clone, Debug)]
pub struct GenericCyclotomic {
    // (exp, (numerator, denominator))
    pub exp_coeffs: HashMap<Z, (i64, u64)>,
    pub order: Z,
}

#[macro_export]
macro_rules! field_axiom_tests {
    ($type: ident) => {
        #[cfg(test)]
        mod tests {
            use super::*;

            #[quickcheck]
            fn zero_is_add_identity(z: $type) -> bool {
                z.clone()
                    .add(&mut $type::zero_order(z.order.clone()))
                    .eq(&mut z.clone())
            }

            #[quickcheck]
            fn add_is_associative(x: $type, y: $type, z: $type) -> bool {
                (x.clone().add(&mut y.clone()))
                    .add(&mut z.clone())
                    .eq(x.clone().add(y.clone().add(&mut z.clone())))
            }

            #[quickcheck]
            fn add_is_commutative(x: $type, y: $type) -> bool {
                x.clone()
                    .add(&mut y.clone())
                    .eq(y.clone().add(&mut x.clone()))
            }

            #[quickcheck]
            fn one_is_mul_identity(z: $type) -> bool {
                let mut same = z
                    .clone()
                    .mul(&mut $type::one_order(z.order.clone()))
                    .clone();
                same.eq(&mut z.clone())
            }

            #[quickcheck]
            fn add_has_inverses(z: $type) -> bool {
                z.clone()
                    .add(z.clone().add_invert())
                    .eq(&mut $type::zero_order(z.order))
            }

            #[quickcheck]
            fn zero_kills_all(z: $type) -> bool {
                $type::zero_order(z.order.clone())
                    .mul(&mut z.clone())
                    .eq(&mut $type::zero_order(z.order))
            }

            #[quickcheck]
            fn mul_is_commutative(x: $type, y: $type) -> bool {
                x.clone()
                    .mul(&mut y.clone())
                    .eq(y.clone().mul(&mut x.clone()))
            }

            #[quickcheck]
            fn mul_is_associative(x: $type, y: $type, z: $type) -> bool {
                (x.clone().mul(&mut y.clone()))
                    .mul(&mut z.clone())
                    .eq(x.clone().mul(y.clone().mul(&mut z.clone())))
            }

            #[quickcheck]
            fn mul_has_inverses(arb_z: $type) -> bool {
                let z = convert_to_base(&arb_z);
                if is_zero(&z) {
                    return true;
                }
                let mut prod = z.clone().mul_invert().mul(&mut z.clone()).clone();
                prod.eq(&mut $type::one_order(z.order))
            }

            #[quickcheck]
            fn mul_distributes_over_add(x: $type, y: $type, z: $type) -> bool {
                x.clone().mul(y.clone().add(&mut z.clone())).eq(x
                    .clone()
                    .mul(&mut y.clone())
                    .add(x.clone().mul(&mut z.clone())))
            }
        }
    };
}

/// Sparse implementation using hash maps.
pub mod sparse;

/// Sparse hash map implementation supporting arbitrarily sized exponents.
pub mod big_sparse;

/// Dense representation using a vector of coefficients.
pub mod dense;

/// Dense representation using structure constants to multiply. Only supports
/// multiplication between elements expressed as elements of the same field,
/// which is the common case.
pub mod structure;

/// Common utilities
pub mod util;

/// Matrix and vector operations.
pub mod linear_algebra;

/// Trait for types that can be used as an exponent. Used mainly for different
/// big integer implementations, and machine-sized integers.
pub mod exponent;