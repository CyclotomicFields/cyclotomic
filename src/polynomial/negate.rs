use std::ops::Neg;
use crate::polynomial::polynomial::Polynomial;

impl Neg for Polynomial {
    type Output = Polynomial;

    fn neg(self) -> Self::Output {
        (&self).neg()
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;
}