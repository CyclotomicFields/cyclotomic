use num::{One, Zero};
use std::ops::{Add, Mul};
use quickcheck::Arbitrary;

type Z = num::bigint::BigInt;
type Q = num::rational::BigRational;
type Exponent = u64;

/// Extremely simple naive implementation. Lots of copying, very inefficient,
/// but algorithmically simple.
mod naive;

/// Provides common field operations. Note that this trait requires both the
/// standard arithmetic traits and versions that allow mutation of the
/// arguments. This is for efficiency reasons, since many algorithms benefit
/// from being able to reduce and expand terms in-place, without changing which
/// element of the field is being represented.
pub trait FieldElement: Eq + PartialEq + Add + Zero + Mul + One {
    fn add_mut(&mut self, other: &mut Self) -> Self;
    fn mul_mut(&mut self, other: &mut Self) -> Self;
    fn inv(&self) -> Self;
    fn inv_mut(&mut self) -> Self;
}

/// Provides convenience functions specific to cyclotomic fields.
pub trait CyclotomicFieldElement: FieldElement {
    /// Returns $\zeta_n$^k.
    fn e(n: u64, k: u64) -> Self;

    /// Multiplies by this scalar. Recall that $\mathbb{Q}(\zeta_n)$ is a
    /// $\mathbb{Q}$-vector space.
    fn scalar_mul(&self, scalar: Q) -> Self;
}

// TODO: work out a way to write tests for a trait then instantiate them instead of
// tests individually for each module.

