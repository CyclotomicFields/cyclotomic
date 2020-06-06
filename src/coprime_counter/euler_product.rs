use crate::coprime_counter::coprime_counter::CoprimeCounter;

type ZPlus = u64;

struct EulerProduct {}

impl EulerProduct {
    fn new() -> EulerProduct {
        EulerProduct {}
    }
}

impl CoprimeCounter for EulerProduct {
    /*
    phi(n) = n * product(1 - (1 / p)) for all p which are the distinct prime
                                      factors of n.

    To compute this then, we need to first get the distinct prime factors of n,
    and then do a product calculation. If n is prime, then of course its only
    prime factor will be itself, and the resulting product will give n - 1.
    */
    fn phi(&self, n: ZPlus) -> ZPlus {
        0
    }
}

#[cfg(test)]
mod phi_tests {
    use super::*;

    #[test]
    fn test_phi() {
        let euler_product = EulerProduct::new();
        assert_eq!(euler_product.phi(1), 0);
    }
}