use std::ops::Mul;

use ::num::integer::lcm;
use num::{One, Zero};

use crate::polynomial::polynomial::{Polynomial, Q, Z};

impl Polynomial {
    pub fn new(coefficients: Vec<Q>) -> Polynomial {
        Polynomial { coefficients }
    }

    pub fn from_small_fractions(coefficient_numerators: Vec<i64>,
                                coefficient_denominators: Vec<i64>) -> Polynomial {
        let mut qs = vec![];
        for i in 0..coefficient_numerators.len() {
            qs.push(Q::new(
                Z::from(coefficient_numerators[i]),
                Z::from(coefficient_denominators[i])));
        }
        Polynomial::new(qs)
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

impl From<Vec<Z>> for Polynomial {
    fn from(coefficients: Vec<Z>) -> Self {
        Polynomial::new(coefficients.iter().map(|c| Q::from(c.clone())).collect())
    }
}

impl From<Vec<i64>> for Polynomial {
    fn from(vec: Vec<i64>) -> Self {
        Polynomial::new(vec.iter().map(|&i| Q::from(Z::from(i))).collect::<Vec<Q>>())
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;
}