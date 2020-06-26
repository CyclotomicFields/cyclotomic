use std::ops::Mul;

use ::num::integer::lcm;
use num::{One, Zero};

use crate::polynomial::polynomial::{Polynomial, Q, Z};

impl Polynomial {
    pub fn new(coefficients: Vec<Z>) -> Polynomial {
        Polynomial { coefficients }
    }

    pub fn truncate_coefficients<T: Zero>(coefficients: &mut Vec<T>) {
        while coefficients.len() > 0 && coefficients[coefficients.len() - 1].is_zero() {
            coefficients.truncate(coefficients.len() - 1);
        }
        if coefficients.is_empty() {
            coefficients.push(T::zero());
        }
    }
}

impl From<Vec<Q>> for Polynomial {
    fn from(coefficients: Vec<Q>) -> Self {
        let mut accumulator = Z::one();
        for c in &coefficients {
            accumulator = lcm(accumulator, c.denom().clone());
        }
        let integer_coefficients = coefficients.iter().map(|c| c.mul(&accumulator).to_integer()).collect();
        Polynomial { coefficients: integer_coefficients }
    }
}

impl From<Vec<i64>> for Polynomial {
    fn from(vec: Vec<i64>) -> Self {
        Polynomial::new(vec.iter().map(|&i| Z::from(i)).collect::<Vec<Z>>())
    }
}

// TODO: Use quickcheck instead
#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_initialise_with_rational_coefficients() {
        #[inline]
        fn vec_q(numerators: Vec<i64>, denominators: Vec<i64>) -> Vec<Q> {
            let mut qs = vec![];
            for i in 0..numerators.len() {
                qs.push(Q::new(Z::from(numerators[i]), Z::from(denominators[i])));
            }
            qs
        }

        assert_eq!(Polynomial::from(vec_q(vec![-7, 2, 1], vec![3, 5, 1])),
                   (Polynomial::from(vec![-35, 6, 15])));

        assert_eq!(Polynomial::from(vec_q(vec![-7, 2, 1, 2, -4], vec![2, 5, 1, -9, 3])),
                   (Polynomial::from(vec![-315, 36, 90, -20, -120])));
    }
}