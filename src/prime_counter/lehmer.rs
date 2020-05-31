extern crate combinations;

use std::ops::{Div, Range};

use crate::prime_counter::prime_counter::PrimeCounter;

use self::combinations::Combinations;
use crate::primes::Primes;
use crate::prime_counter::legendre::Legendre;

type R = f64;
type ZPlus = u64;
type Z = i64;

pub struct Lehmer<'a> {
    primes: &'a Primes,
    legendre: Legendre<'a>,
}

impl<'a> PrimeCounter for Lehmer<'a> {
    fn pi(&self, x: R) -> ZPlus {
        self.pi_prime(x, &self.primes)
    }
}

impl<'a> Lehmer<'a> {
    pub fn new(primes: &Primes) -> Lehmer {
        Lehmer { primes, legendre: Legendre::new(primes) }
    }

    fn pi_prime(&self, x: R, relevant_primes: &Primes) -> ZPlus {
        return if x < 2.0 {
            0
        } else if x < 3.0 {
            1
        } else if x.floor() == 3.0 {
            2
        } else {
            /*
            Lehmer's method is an extension of Legendre's method. It scales
            better for larger values of x, because it takes a Legendre sum that
            only goes up to pi(x ^ 1/4), rather than for Legendre's method,
            which goes up to pi(sqrt(x)). To account for this, it adds a couple
            more terms.
            */
            let fourth_root_x = x.powf(1.0 / 4.0).floor();
            let a = relevant_primes.pi(fourth_root_x);
            let cube_root_x = x.powf(1.0 / 3.0).floor();
            let c = relevant_primes.pi(cube_root_x);
            let square_root_x = x.sqrt().floor();
            let b = relevant_primes.pi(square_root_x);

            let square_root_primes = relevant_primes.range(0, square_root_x as ZPlus).unwrap();
            if let Some(fourth_root_primes) = relevant_primes.range(0, fourth_root_x as ZPlus) {
                let legendre_sum = self.legendre.legendre_sum(x, fourth_root_primes.to_vec());
                let constant_term = ((((b + a) as R - 2.0) * ((b - a) as R + 1.0)) / 2.0).floor() as Z;
                let one_prime_pi_sum = relevant_primes.range(fourth_root_x as ZPlus, square_root_x as ZPlus).unwrap().to_vec()
                    .iter()
                    .map(|&p| self.pi_prime(x / p as R, &square_root_primes))
                    .sum::<ZPlus>();
                let two_prime_pi_sum = ((a + 1)..(c + 1))
                    .map(|i| (i..(1 + relevant_primes.pi((x / relevant_primes.nth(i as usize).unwrap() as R).sqrt())))
                        .map(|j| {
                            let p_i = relevant_primes.nth(i as usize).unwrap();
                            let p_j = relevant_primes.nth(j as usize).unwrap();
                            self.pi_prime(x / (p_i * p_j) as R, &square_root_primes) - (j as ZPlus - 1)
                        })
                        .sum::<ZPlus>())
                    .sum::<ZPlus>();

                return (legendre_sum + constant_term - one_prime_pi_sum as Z - two_prime_pi_sum as Z) as ZPlus;
            } else {
                panic!("Not enough in-memory primes to compute pi({})", x);
            }
        };
    }
}

#[cfg(test)]
mod lehmer_tests {
    use std::path::Path;

    use crate::prime_table::PrimeTableReader;

    use super::*;

    #[test]
    fn test_lehmer_fast() {
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199]);
        let strategy: Lehmer = Lehmer::new(&primes);
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
        assert_eq!(strategy.pi(10000.0_f64), 1229);
    }

    #[test]
    fn test_lehmer_single() {
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]);
        let strategy: Lehmer = Lehmer::new(&primes);
        assert_eq!(strategy.pi(200.0_f64), 46);
    }

    #[test]
    fn test_lehmer_slow() {
        if let Some(prime_table_reader) = PrimeTableReader::first_million_from_file() {
            let primes = Primes::new(prime_table_reader.first_million_primes());
            let strategy: Lehmer = Lehmer::new(&primes);
            assert_eq!(strategy.pi(20000000.0_f64), 1270607);
        }
    }
}