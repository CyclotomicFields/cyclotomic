extern crate combinations;

use crate::prime_counter::prime_counter::PrimeCounter;

use std::ops::Div;

use self::combinations::Combinations;

type R = f64;
type ZPlus = u64;
type Z = i64;

pub struct Legendre {
    primes: Vec<ZPlus>
}

impl PrimeCounter for Legendre {
    fn pi(&self, x: R) -> ZPlus {
        self.pi_prime(x, self.primes.clone())
    }
}

impl Legendre {
    pub fn new(primes: Vec<ZPlus>) -> Legendre {
        Legendre { primes }
    }

    /*
    Legendre's method in a sentence:

    The number of primes below x equals the total number of numbers below x minus the number of
    composites below x.

    The formula looks like this:

    pi(x) = -1 + x + pi(sqrt(x))
            - (floor each (x div each-right primes below sqrt(x))
            + (floor each (x div each-right products of all pairs of primes below sqrt(x))
            - (floor each (x div each-right products of all groups of three primes below sqrt(x))
            + (floor each (x div each-right products of all groups of four primes below sqrt(x))
            ...
    */
    fn pi_prime(&self, x: R, relevant_primes: Vec<ZPlus>) -> ZPlus {
        return if x < 2.0 {
            0
        } else if x == 2.0 {
            1
        } else if x == 3.0 {
            2
        } else {
            let mut relevant_primes: Vec<ZPlus> = relevant_primes.clone();
            relevant_primes.retain(|&p| p <= x.sqrt().floor() as ZPlus);
            println!("Finding pi({}) using primes: {:?}", x, relevant_primes);

            /*
            This inline function will divide x by the product of the given prime numbers, to
            find a term for the number of composite numbers that correspond to this group of
            primes. The reasoning for this is given just above the loop that sums up these terms.
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
            The terms for the composite number component have different signs. This is because of
            the way the formula accounts for overcompensation in earlier terms in its expansion.
            More details are in the comment above the loop that adds up these terms. That comment
            goes into more detail on the composite number count terms.
            */
            #[inline]
            fn composite_term_sign(number_of_primes_in_term_denominator: usize) -> Z {
                return if number_of_primes_in_term_denominator % 2 == 0 { -1 } else { 1 };
            }

            /*
            All composite numbers in the interval [1, x] have at least one prime factor less than
            or equal to sqrt(x), so if we take the sum of all values of x/p for these primes that
            should give it to us. In doing this, we don't want to consider the terms 1*p for
            prime p as composites. This is why we subtract the term pi(sqrt(x)).

            Of course, some of the composites in [1, x] are divisible by two primes!
            These composites will correspondingly have been accounted for twice, in the
            previous calculation. We will need to add another corrective factor in the other
            direction for these composite numbers. Of course this correction then causes an
            inaccuracy regarding those composites with three prime factors, and so on.

            This is why we have all these corrective terms. It is also why the terms based on
            even-numbered groups of primes have one sign, and odd-numbered groups of primes have the
            other.
            */

            let x_z = x as Z;
            let mut composite_number_count = 0;
            for i in 1..relevant_primes.len() {
                let sign = composite_term_sign(i);
                composite_number_count += sign * Combinations::new(relevant_primes.clone(), i)
                    .map(|prime_product_group| composite_numbers_term(x_z, &prime_product_group))
                    .sum::<Z>();
            }
            let sign = composite_term_sign(relevant_primes.len());
            composite_number_count += sign * composite_numbers_term(x_z, &relevant_primes);
            composite_number_count -= self.pi_prime(x.sqrt().floor(), relevant_primes) as Z;

            /*
            Legendre's method is based on the fact that the number of all integers below x
            (= -1 + x) equals the number of primes (= pi(x), the result) plus the number of
            composite numbers. The terms accounting for the composite numbers are the most complex,
            which is detailed above.
            */
            return ((-1 + x as Z) - composite_number_count) as ZPlus;
        };
    }
}

#[cfg(test)]
mod legendre_tests {
    use std::path::Path;

    use crate::prime_table::PrimeTableReader;

    use super::*;

    #[test]
    fn test_legendre_fast() {
        let strategy: Legendre = Legendre::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]);
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
        let strategy: Legendre = Legendre::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]);
        assert_eq!(strategy.pi(5000.0_f64), 669);
    }
}