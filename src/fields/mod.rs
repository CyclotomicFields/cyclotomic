use num::{One, Zero};
use quickcheck::Arbitrary;

type Z = num::BigInt;
type Q = num::rational::BigRational;

/// Sparse implementation using hash maps.
pub mod sparse;

pub trait AdditiveGroup {
    fn add_invert(&mut self) -> &mut Self;

    /// Adds z to self in place, so self = self + z
    fn add(&mut self, z: &mut Self) -> &mut Self;
}

pub trait MultiplicativeGroup {
    /// Inverts self in place.
    fn mul_invert(&mut self) -> &mut Self;

    /// Multiplies self by z in place, so self = self * z
    fn mul(&mut self, z: &mut Self) -> &mut Self;
}

/// Provides operations for fields. Expected to satisfy the field axioms.
pub trait FieldElement: AdditiveGroup + MultiplicativeGroup {
    /// Equality, but can shuffle coefficients around and simplify expressions
    /// for greater efficiency.
    fn eq(&mut self, other: &mut Self) -> bool;
}

/// Provides convenience functions specific to cyclotomic fields.
pub trait CyclotomicFieldElement: FieldElement {
    /// Returns $\zeta_n$^k.
    fn e(n: i64, k: i64) -> Self;

    /// Multiplies in-place by scalar. Recall that $\mathbb{Q}(\zeta_n)$ is a
    /// $\mathbb{Q}$-vector space.
    fn scalar_mul(&mut self, scalar: &Q) -> &mut Self;
}

// TODO: work out a way to write tests for a trait then instantiate them instead of
// tests individually for each module.

