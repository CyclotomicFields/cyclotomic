use crate::coprime_counter::coprime_counter::CoprimeCounter;
use crate::prime_factors::prime_factorize::PrimeFactorize;
use crate::divisors::library_divisors::LibraryDivisors;
use crate::prime_factors::recursive_prime_factorize::RecursivePrimeFactorize;
use num::ToPrimitive;

type ZPlus = u64;
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
    fn phi(&self, n: ZPlus) -> ZPlus {
        return if n == 0 {
            panic!("Can't evaluate phi(0)")
        } else if n == 1 {
            1
        } else if n % 2 == 0 {
            let n_over_two = (n / 2);
            if n_over_two % 2 == 0 {
                2 * self.phi(n_over_two)
            } else {
                self.phi(n_over_two)
            }
        } else {
            let mut prime_factors = self.prime_factorizer.prime_factors(Z::from(n));
            prime_factors.dedup();
            (n as R * prime_factors
                .iter()
                .map(|p| 1.0 - (p.to_f64().unwrap()).recip())
                .product::<R>()) as ZPlus
        };
    }
}

#[cfg(test)]
mod phi_tests {
    use super::*;

    #[test]
    fn test_phi_small() {
        let euler_product = HybridCoprimeCounter::default();
        assert_eq!(euler_product.phi(1), 1);
        assert_eq!(euler_product.phi(2), 1);
        assert_eq!(euler_product.phi(8), 4);
        assert_eq!(euler_product.phi(3), 2);
        assert_eq!(euler_product.phi(9), 6);
        assert_eq!(euler_product.phi(100), 40);
        assert_eq!(euler_product.phi(370416), 123456);
    }

    #[test]
    fn test_phi_medium() {
        let euler_product = HybridCoprimeCounter::default();
        // It's not actually accurate for numbers much bigger than this.
        assert_eq!(euler_product.phi(123456789101112), 41135376570624);
        assert_eq!(euler_product.phi(1234567891011121), 1126818037342720);
        assert_eq!(euler_product.phi(12345678910111212), 4055421449693376);
        assert_eq!(euler_product.phi(12345678910111), 12345678910110);
    }
}