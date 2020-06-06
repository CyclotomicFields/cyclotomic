use crate::coprime_counter::coprime_counter::CoprimeCounter;
use crate::prime_factors::prime_factorize::PrimeFactorize;
use crate::divisors::library_divisors::LibraryDivisors;
use crate::prime_factors::recursive_prime_factorize::RecursivePrimeFactorize;
use num::{ToPrimitive, Zero, One, Integer};
use std::ops::Div;

type R = f64;
type Q = num::rational::BigRational;
type Z = num::bigint::BigInt;

struct HybridCoprimeCounter<PF: PrimeFactorize> {
    prime_factorizer: PF
}

impl<PF: PrimeFactorize> HybridCoprimeCounter<PF> {
    pub fn new(prime_factorizer: PF) -> HybridCoprimeCounter<PF> {
        HybridCoprimeCounter { prime_factorizer }
    }
}

impl HybridCoprimeCounter<RecursivePrimeFactorize<LibraryDivisors>> {
    pub fn default() -> HybridCoprimeCounter<RecursivePrimeFactorize<LibraryDivisors>> {
        HybridCoprimeCounter::new(RecursivePrimeFactorize::default())
    }
}

impl<PF: PrimeFactorize> CoprimeCounter for HybridCoprimeCounter<PF> {
    /*
    The base case of the implementation: Euler's product formula

    phi(n) = n * product(1 - (1 / p)) for all p which are the distinct prime
                                      factors of n.

    To compute this then, we need to first get the distinct prime factors of n,
    and then do a product calculation. If n is prime, then of course its only
    prime factor will be itself, and the resulting product will give n - 1.


    Reduction into divisors

    phi(m * n) = phi(m) * phi(n) * d / phi(d) where d = gcd(m, n)


    Special case: Reduction by halving

    phi(2m) = 2 * phi(m) if m is even
    phi(2m) = phi(m)     if m is odd

    So we can always reduce the size of n by halving it until we get an odd
    integer.

    */
    fn phi(&self, n: &Z) -> Z {
        return if n.is_zero() {
            panic!("Can't evaluate phi(0)")
        } else if n.is_one() {
            Z::one()
        } else if n.is_even() {
            let n_over_two = (n.div(&Z::from(2)));
            if n_over_two.is_even() {
                2 * self.phi(&n_over_two)
            } else {
                self.phi(&n_over_two)
            }
        } else {
            let mut prime_factors = self.prime_factorizer.prime_factors(n.clone());
            prime_factors.dedup();
            Z::from((n.to_f64().unwrap() * prime_factors
                .iter()
                .map(|p| 1.0 - (p.to_f64().unwrap()).recip())
                .product::<R>()) as i64)
        };
    }
}

#[cfg(test)]
mod phi_tests {
    use super::*;

    #[test]
    fn test_phi_small() {
        let euler_product = HybridCoprimeCounter::default();
        assert_eq!(euler_product.phi(&Z::from(1)), Z::from(1));
        assert_eq!(euler_product.phi(&Z::from(2)), Z::from(1));
        assert_eq!(euler_product.phi(&Z::from(8)), Z::from(4));
        assert_eq!(euler_product.phi(&Z::from(3)), Z::from(2));
        assert_eq!(euler_product.phi(&Z::from(9)), Z::from(6));
        assert_eq!(euler_product.phi(&Z::from(100)), Z::from(40));
        assert_eq!(euler_product.phi(&Z::from(370416)), Z::from(123456));
    }

    #[test]
    fn test_phi_medium() {
        let euler_product = HybridCoprimeCounter::default();
        // It's not actually accurate for numbers much bigger than this.
        assert_eq!(euler_product.phi(&Z::from(123456789101112_i64)), Z::from(41135376570624_i64));
        assert_eq!(euler_product.phi(&Z::from(1234567891011121_i64)), Z::from(1126818037342720_i64));
        assert_eq!(euler_product.phi(&Z::from(12345678910111212_i64)), Z::from(4055421449693376_i64));
        assert_eq!(euler_product.phi(&Z::from(12345678910111_i64)), Z::from(12345678910110_i64));
    }
}