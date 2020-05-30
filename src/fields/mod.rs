use std::ops::{Add, Mul};
use num::{Zero, One};

/// Extremely simple naive implementation. Lots of copying, very inefficient
/// but algorithmically simple.
mod naive;

/// Provides common field operations. Note that this trait requires both the
/// standard arithmetic traits, and versions that are add to an existing
/// element in-place. Note that, for efficiency reasons, the functions are
/// allowed to mutate both arguments, but should not change the field element
/// represented by `other`. This is not enforced at the type level.
pub trait FieldElement: Eq + PartialEq + Add + Zero + Mul + One {
    fn add_mut(&mut self, other: &mut Self);
    fn mul_mut(&mut self, other: &mut Self);
    fn inv(&self) -> Self;
    fn inv_mut(&mut self);
}

