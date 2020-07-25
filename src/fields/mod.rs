use num::{One, Zero};
use quickcheck::Arbitrary;

pub type Z = num::BigInt;
pub type Q = num::rational::BigRational;

pub trait AdditiveGroupElement {
    /// Adds z to self in place, so self = self + z
    fn add(&mut self, z: &mut Self) -> &mut Self;

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
pub trait CyclotomicFieldElement: FieldElement + Clone {
    /// Returns $\zeta_n$^k.
    fn e(n: i64, k: i64) -> Self;

    /// Multiplies in-place by scalar. Recall that $\mathbb{Q}(\zeta_n)$ is a
    /// $\mathbb{Q}$-vector space.
    fn scalar_mul(&mut self, scalar: &Q) -> &mut Self;

    /// Gives zero expressed as an element of $\mathbb{Q}(\zeta_n)$
    fn zero_order(n: i64) -> Self;

    /// Gives one expressed as an element of $\mathbb{Q}(\zeta_n)$
    fn one_order(n: i64) -> Self;
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

/// Dense representation using a vector of coefficients.
pub mod dense;

