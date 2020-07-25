use crate::polynomial::polynomial::Polynomial;
use std::ops::Neg;

impl Polynomial {
    pub fn neg(&self) -> Polynomial {
        Polynomial::new(self.coefficients.iter().map(|c| c.neg()).collect())
    }
}

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
