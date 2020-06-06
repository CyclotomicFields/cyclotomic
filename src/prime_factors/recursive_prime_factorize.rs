extern crate num;

use self::num::{BigInt, BigRational, Zero, ToPrimitive};

use crate::prime_factors::prime_factorize::PrimeFactorize;
use crate::divisors::divisors::Divisors;
use crate::divisors::library_divisors::LibraryDivisors;

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
    fn prime_factors(&self, n: u64) -> Vec<u64> {
        let divisors = self.divisors_strategy.divisors_without_one(Z::from(n));
        return if divisors.len() == 1 {
            let d = divisors[0].clone();
            if let Some(d_zplus) = d.to_u64() {
                vec![d_zplus]
            } else {
                panic!("Couldn't convert {} into ZPlus", d);
            }
        } else {
            let d = divisors[0].clone();
            if let Some(d_zplus) = d.to_u64() {
                let mut result = vec![d_zplus];
                result.extend(self.prime_factors(n / d_zplus));
                return result;
            } else {
                panic!("Couldn't convert {} into ZPlus", d);
            }
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prime_factorize() {
        let factorizer = RecursivePrimeFactorize::default();
        assert_eq!(factorizer.prime_factors(2), vec![2]);
        assert_eq!(factorizer.prime_factors(3), vec![3]);
        assert_eq!(factorizer.prime_factors(4), vec![2, 2]);
        assert_eq!(factorizer.prime_factors(12), vec![2, 2, 3]);
        assert_eq!(factorizer.prime_factors(130), vec![2, 5, 13]);
        assert_eq!(factorizer.prime_factors(1300), vec![2, 2, 5, 5, 13]);
        assert_eq!(factorizer.prime_factors(67), vec![67]);
    }
}