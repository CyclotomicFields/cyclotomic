extern crate num;

use std::ops::Div;
use num::BigInt;
use crate::divisors::divisors::Divisors;
use crate::divisors::library_divisors::LibraryDivisors;
use crate::prime_factors::prime_factorize::PrimeFactorize;

use self::num::{One};

type Z = num::bigint::BigInt;

pub struct RecursivePrimeFactorize<D: Divisors> {
    divisors_strategy: D
}

impl<D: Divisors> RecursivePrimeFactorize<D> {
    pub fn new(divisors_strategy: D) -> RecursivePrimeFactorize<D> {
        RecursivePrimeFactorize { divisors_strategy }
    }
}

impl RecursivePrimeFactorize<LibraryDivisors> {
    pub fn default() -> RecursivePrimeFactorize<LibraryDivisors> {
        RecursivePrimeFactorize::new(LibraryDivisors::new())
    }
}

impl<D: Divisors> PrimeFactorize for RecursivePrimeFactorize<D> {
    fn prime_factors(&self, n: &Z) -> Vec<Z> {
        if n.is_one() {
            return vec![];
        }
        let divisors = self.divisors_strategy.divisors_without_one(n);
        let d = divisors[0].clone();
        let mut result = vec![d.clone()];
        return if divisors.len() == 1 {
            result
        } else {
            result.extend(self.prime_factors(&n.div(&d)));
            result
        };
    }
}

#[cfg(test)]
mod recursive_prime_factorize_tests {
    use super::*;

    #[test]
    fn test_prime_factorize() {
        let factorizer = RecursivePrimeFactorize::default();
        assert_eq!(factorizer.prime_factors(&Z::from(1)), vec_z(vec![]));
        assert_eq!(factorizer.prime_factors(&Z::from(2)), vec_z(vec![2]));
        assert_eq!(factorizer.prime_factors(&Z::from(3)), vec_z(vec![3]));
        assert_eq!(factorizer.prime_factors(&Z::from(4)), vec_z(vec![2, 2]));
        assert_eq!(factorizer.prime_factors(&Z::from(12)), vec_z(vec![2, 2, 3]));
        assert_eq!(factorizer.prime_factors(&Z::from(130)), vec_z(vec![2, 5, 13]));
        assert_eq!(factorizer.prime_factors(&Z::from(1300)), vec_z(vec![2, 2, 5, 5, 13]));
        assert_eq!(factorizer.prime_factors(&Z::from(67)), vec_z(vec![67]));
    }

    fn vec_z(vec: Vec<i64>) -> Vec<Z> {
        vec.iter().map(|&i| BigInt::from(i)).collect()
    }
}