use num::One;

pub type Z = num::bigint::BigInt;
pub type Q = num::rational::BigRational;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coefficients: Vec<Q>,
}

impl Polynomial {
    pub fn degree(&self) -> usize {
        return self.coefficients.len() - 1;
    }

    pub fn leading_term_coefficient(&self) -> &Q {
        &self.coefficients[self.coefficients.len() - 1]
    }

    pub fn constant_term_coefficient(&self) -> &Q {
        &self.coefficients[0]
    }

    pub fn is_monic(&self) -> bool {
        return self.leading_term_coefficient().is_one();
    }
}

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        self.coefficients == other.coefficients
    }
}

// TODO: Use quickcheck instead
#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_monic_check() {
        // p = 1
        assert!(Polynomial::from(vec![1]).is_monic());
        // p = 2t + 1
        assert!(!Polynomial::from(vec![1, 2]).is_monic());
        // p = t^5 - t^3 - 2t^2 - 1
        assert!(Polynomial::from(vec![-1, 0, -2, -1, 0, 1]).is_monic());
    }
}
