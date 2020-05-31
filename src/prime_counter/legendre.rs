extern crate combinations;

use std::ops::Div;

use crate::prime_counter::prime_counter::PrimeCounter;

use self::combinations::Combinations;
use crate::primes::Primes;

type R = f64;
type ZPlus = u64;
type Z = i64;

pub struct Legendre<'a> {
    primes: &'a Primes
}

impl<'a> PrimeCounter for Legendre<'a> {
    fn pi(&self, x: R) -> ZPlus {
        self.pi_prime(x, &self.primes)
    }
}

impl<'a> Legendre<'a> {
    pub fn new(primes: &Primes) -> Legendre {
        Legendre { primes }
    }

    /*
    Legendre's method in a sentence:

    The number of primes below x equals the total number of numbers below x
    minus the number of composites below x.

    The formula looks like this:

    pi(x) = -1 + x + pi(sqrt(x))
            - (floor each (x div each-right primes below sqrt(x))
            + (floor each (x div each-right products of all pairs of primes below sqrt(x))
            - (floor each (x div each-right products of all groups of three primes below sqrt(x))
            + (floor each (x div each-right products of all groups of four primes below sqrt(x))
            ...
    */
    fn pi_prime(&self, x: R, relevant_primes: &Primes) -> ZPlus {
        return if x < 2.0 {
            0
        } else if x < 3.0 {
            1
        } else if x.floor() == 3.0 {
            2
        } else {
            if let Some(primes) = relevant_primes.range(0, x.sqrt().floor() as ZPlus) {
                let primes_vec = primes.to_vec();
                /*
                Legendre's method is based on the fact that the number of all
                integers below x (= -1 + x) equals the number of primes
                (= pi(x), the result) plus the number of composite numbers.
                */
                return (self.legendre_sum(x, primes_vec) + self.pi_prime(x.sqrt().floor(), &primes) as Z - 1) as ZPlus;
            } else {
                panic!("Not enough in-memory primes to compute pi({})", x);
            }
        };
    }

    /*
    This sum is used for the Legendre, Meissel and Lehmer methods for
    calculating pi(x). The sum counts the positive integers less than or equal
    to x, not divisible by any one of the primes in the relevant_primes vector
    reference.

    For input x and a set of primes P, the legendre sum is equal to:

    floor(x) - sum[over P](x / p_i) + sum[over P](x / p_i*p_j) - sum[over P](x / p_i*p_j*p_k) ...

    It is denoted with a capital Phi as a function Phi(x, n), which yields the
    legendre sum of x using the first n primes.
    */
    pub fn legendre_sum(&self, x: R, relevant_primes: &Vec<ZPlus>) -> Z {
        /*
        If there are no relevant primes, then the result will be equal to floor(x)
        */
        if relevant_primes.is_empty() {
            return x.floor() as Z;
        }

        /*
        This inline function will divide x by the product of the given prime
        numbers, to find a term for the number of composite numbers that
        correspond to this group of primes. The reasoning for this is given
        just above the loop that sums up these terms.
        */
        #[inline]
        fn composite_numbers_term(x: Z, prime_product_group: &Vec<ZPlus>) -> Z {
            let mut prime_product = 1;
            for &q in prime_product_group {
                prime_product *= q as u128;
            }
            return x.div(prime_product as Z);
        }

        /*
        The terms for the composite number component have different signs. This
        is because of the way the formula accounts for overcompensation in
        earlier terms in its expansion. More details are in the comment above
        the loop that adds up these terms. That comment goes into more detail
        on the composite number count terms.
        */
        #[inline]
        fn composite_term_sign(number_of_primes_in_term_denominator: usize) -> Z {
            return if number_of_primes_in_term_denominator % 2 == 0 { -1 } else { 1 };
        }

        /*
        All composite numbers in the interval [1, x] have at least one prime
        factor less than or equal to sqrt(x), so if we take the sum of all
        values of x/p for these primes that should give it to us. In doing
        this, we don't want to consider the terms 1*p for prime p as composites.
        This is why we subtract the term pi(sqrt(x)).

        Of course, some of the composites in [1, x] are divisible by two primes!
        These composites will correspondingly have been accounted for twice, in
        the previous calculation. We will need to add another corrective factor
        in the other direction for these composite numbers. Of course this
        correction then causes an inaccuracy regarding those composites with
        three prime factors, and so on.

        This is why we have all these corrective terms. It is also why the terms
        based on even-numbered groups of primes have one sign, and odd-numbered
        groups of primes have the other.
        */
        let x_z = x as Z;
        let mut legendre_sum = x_z;
        for i in 1..relevant_primes.len() {
            let sign = composite_term_sign(i);
            legendre_sum -= sign * Combinations::new(relevant_primes.clone(), i)
                .map(|prime_product_group| composite_numbers_term(x_z, &prime_product_group))
                .sum::<Z>();
        }
        let sign = composite_term_sign(relevant_primes.len());
        legendre_sum -= sign * composite_numbers_term(x_z, &relevant_primes);
        return legendre_sum;
    }

    /*
    This method is for convenience, to match the mathematical notation.

    The result is the number of positive integers less than or equal
    to x, not divisible by any of the first a primes.
    */
    pub fn legendre_sum_phi(&self, x: R, a: ZPlus) -> Z {
        if let Some(primes) = self.primes.first_n(a) {
            return self.legendre_sum(x, primes.to_vec());
        } else {
            panic!("Couldn't evaluate phi({}, {}) using {} primes", x, a, self.primes.len())
        }
    }
}

#[cfg(test)]
mod legendre_tests {
    use std::path::Path;

    use crate::prime_table::PrimeTableReader;

    use super::*;

    #[test]
    fn test_legendre_fast() {
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67]);
        let strategy: Legendre = Legendre::new(&primes);
        assert_eq!(strategy.pi(1.0_f64), 0);
        assert_eq!(strategy.pi(2.0_f64), 1);
        assert_eq!(strategy.pi(3.0_f64), 2);
        assert_eq!(strategy.pi(4.0_f64), 2);
        assert_eq!(strategy.pi(5.0_f64), 3);
        assert_eq!(strategy.pi(10.0_f64), 4);
        assert_eq!(strategy.pi(100.0_f64), 25);
        assert_eq!(strategy.pi(200.0_f64), 46);
        assert_eq!(strategy.pi(1000.0_f64.sqrt()), 11);
        assert_eq!(strategy.pi(1000.0_f64), 168);
        assert_eq!(strategy.pi(2000.0_f64), 303);
        assert_eq!(strategy.pi(3000.0_f64), 430);
        assert_eq!(strategy.pi(4000.0_f64), 550);
    }

    #[test]
    fn test_legendre_slow() {
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]);
        let strategy: Legendre = Legendre::new(&primes);
        assert_eq!(strategy.pi(5000.0_f64), 669);
    }

    #[test]
    fn test_legendre_sum() {
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]);
        let strategy: Legendre = Legendre::new(&primes);
        let no_primes = vec![];
        let seven_primes = vec![2, 3, 5, 7, 11, 13, 17];
        assert_eq!(strategy.legendre_sum(4.0, &no_primes), 4);
        assert_eq!(strategy.legendre_sum(7.0, &no_primes), 7);
        assert_eq!(strategy.legendre_sum(15.6, &no_primes), 15);
        assert_eq!(strategy.legendre_sum(351.12452, &no_primes), 351);
        assert_eq!(strategy.legendre_sum(100000.0, &seven_primes), 18053);
    }

    #[test]
    fn test_legendre_sum_phi() {
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101]);
        let strategy: Legendre = Legendre::new(&primes);
        assert_eq!(strategy.legendre_sum_phi(100000.0, 7), 18053);
    }
}