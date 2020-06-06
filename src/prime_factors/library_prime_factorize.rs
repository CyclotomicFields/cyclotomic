extern crate num;

use self::num::{BigInt, BigRational, Zero, ToPrimitive};

use crate::prime_factors::prime_factorize::PrimeFactorize;
use crate::divisors::divisors::Divisors;

type Z = num::bigint::BigInt;

// Very lazy but it'll do the job for now!
struct LibraryPrimeFactorize<D: Divisors> {
    divisors_strategy: D
}

impl<D: Divisors> LibraryPrimeFactorize<D> {
    pub fn new(divisors_strategy: D) -> LibraryPrimeFactorize<D> {
        LibraryPrimeFactorize { divisors_strategy }
    }
}

impl<D: Divisors> PrimeFactorize for LibraryPrimeFactorize<D> {
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
    use crate::divisors::library_divisors::LibraryDivisors;

    #[test]
    fn test_prime_factorize() {
        let factorizer = LibraryPrimeFactorize::new(LibraryDivisors::new());
        assert_eq!(factorizer.prime_factors(2), vec![2]);
        assert_eq!(factorizer.prime_factors(3), vec![3]);
        assert_eq!(factorizer.prime_factors(4), vec![2, 2]);
        assert_eq!(factorizer.prime_factors(12), vec![2, 2, 3]);
        assert_eq!(factorizer.prime_factors(130), vec![2, 5, 13]);
    }
}