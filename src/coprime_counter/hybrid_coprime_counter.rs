use crate::coprime_counter::coprime_counter::CoprimeCounter;

type ZPlus = u64;

struct HybridCoprimeCounter {}

impl HybridCoprimeCounter {
    fn new() -> HybridCoprimeCounter {
        HybridCoprimeCounter {}
    }
}

impl CoprimeCounter for HybridCoprimeCounter {
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
        if n % 2 == 0 {
            let n_over_two = (n / 2);
            return if n_over_two % 2 == 0 {
                2 * self.phi(n_over_two)
            } else {
                self.phi(n_over_two)
            };
        }
        1
    }
}

#[cfg(test)]
mod phi_tests {
    use super::*;

    #[test]
    fn test_phi() {
        let euler_product = HybridCoprimeCounter::new();
        assert_eq!(euler_product.phi(1), 1);
        assert_eq!(euler_product.phi(2), 1);
        assert_eq!(euler_product.phi(8), 4);
    }
}