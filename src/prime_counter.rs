extern crate combinations;

use self::combinations::Combinations;
use std::ops::Div;
use std::panic::resume_unwind;

type R = f64;
type Z = i64;
type ZPlus = u64;

pub trait PrimeCounter {
    fn pi(&self, x: R) -> Z;
}

pub struct Legendre {
    primes: Vec<ZPlus>
}

impl Legendre {
    fn new(primes: Vec<ZPlus>) -> Legendre {
        Legendre { primes }
    }
}

impl PrimeCounter for Legendre {
    // pi(x) = -1 + pi(sqrt(x)) + x
    //         - (floor each (x div each-right primes below sqrt(x))
    //         + (floor each (x div each-right products of all pairs of primes below sqrt(x))
    //         - (floor each (x div each-right products of all groups of three primes below sqrt(x))
    //         + (floor each (x div each-right products of all groups of four primes below sqrt(x))
    //         ...
    fn pi(&self, x: R) -> Z {
        return if x < 2.0 {
            0
        } else if x == 2.0 {
            1
        } else {
            println!("\nFinding pi({})", x);
            let sqrt_x = x.sqrt().floor();
            let pi_sqrt_x: ZPlus = self.pi(sqrt_x) as ZPlus;
            println!("sqrt(x) = {}", sqrt_x);
            println!("pi(sqrt(x)) = {}", pi_sqrt_x);
            let mut relevant_primes = self.primes.clone();
            relevant_primes.retain(|&p| p <= sqrt_x as ZPlus);
            println!("Using primes {:?}", relevant_primes);

            let mut result: Z = 0;

            result += pi_sqrt_x as Z;
            result += x as Z;
            result -= 1;

            println!("The result before the corrective terms is {}", result);

            if relevant_primes.is_empty() {
                return result as Z;
            }

            let mut corrective_term_sum = 0;

            fn add_corrective_term(x: R, prime_product_group: Vec<ZPlus>) -> Z {
                return if prime_product_group.len() % 2 == 0 {
                    let mut prime_product = 1;
                    for q in prime_product_group.clone() {
                        prime_product *= q;
                    }
                    x.div(prime_product as R)
                } else {
                    let mut prime_product = 1;
                    for q in prime_product_group.clone() {
                        prime_product *= q;
                    }
                    -x.div(prime_product as R)
                } as Z;
            }

            for i in 1..relevant_primes.len() + 1 {
                println!("Finding prime groups of size {}", i);
                if i == relevant_primes.len() {
                    corrective_term_sum += add_corrective_term(x, relevant_primes.clone());
                } else {
                    for prime_product_group in Combinations::new(relevant_primes.clone(), i) {
                        corrective_term_sum += add_corrective_term(x, prime_product_group);
                    }
                }
            }

            println!("The corrective term is {}", corrective_term_sum);

            return (result + corrective_term_sum) as Z;
        };
    }
}

#[cfg(test)]
mod pi_tests {
    use super::*;

    #[test]
    fn test_combinations() {
        let computed: Vec<_> = Combinations::new(vec![1, 2, 2, 3, 4], 3).collect();
        let expected = vec![
            vec![1, 2, 2],
            vec![1, 2, 3],
            vec![1, 2, 4],
            vec![1, 3, 4],
            vec![2, 2, 3],
            vec![2, 2, 4],
            vec![2, 3, 4],
        ];
        assert!(computed == expected)
    }

    #[test]
    fn test_legendre() {
        let strategy: Legendre = Legendre::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]);
        assert_eq!(strategy.pi(1.0_f64), 0);
        assert_eq!(strategy.pi(2.0_f64), 1);
        assert_eq!(strategy.pi(3.0_f64), 2);
        assert_eq!(strategy.pi(4.0_f64), 2);
        assert_eq!(strategy.pi(5.0_f64), 3);
        assert_eq!(strategy.pi(10.0_f64), 4);
        assert_eq!(strategy.pi(100.0_f64), 25);
        assert_eq!(strategy.pi(200.0_f64), 46);
        assert_eq!(strategy.pi(1000.0_f64.sqrt()), 11);
    }
}