extern crate combinations;

use std::ops::Div;

use crate::prime_counter::prime_counter::PrimeCounter;

use self::combinations::Combinations;
use crate::primes::Primes;
use crate::prime_counter::legendre::Legendre;

type R = f64;
type ZPlus = u64;
type Z = i64;

pub struct Meissel<'a> {
    primes: &'a Primes,
    legendre: Legendre<'a>,
}

impl<'a> PrimeCounter for Meissel<'a> {
    fn pi(&self, x: R) -> ZPlus {
        self.pi_prime(x, &self.primes)
    }
}

impl<'a> Meissel<'a> {
    pub fn new(primes: &Primes) -> Meissel {
        Meissel { primes, legendre: Legendre::new(primes) }
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
            Meissel's method is an extension of Legendre's method. It scales
            better for larger values of x, because it takes a Legendre sum that
            only goes up to pi(cbrt(x)), rather than for Legendre's method,
            which goes up to pi(sqrt(x)). To account for this, it adds a couple
            of other terms - a constant term, and a summation over the values of
            pi(x / p) for the primes in the interval [ cbrt(x), sqrt(x) ).
            */
            let cbrt_x = x.powf(1.0 / 3.0).floor() as ZPlus;
            if let (Some(cbrt_primes), Some(intermediate_primes)) = (relevant_primes.range(0, cbrt_x), relevant_primes.range(cbrt_x, x.sqrt().floor() as ZPlus)) {
                let legendre_sum = self.legendre.legendre_sum(x, cbrt_primes.to_vec());
                let pi_cbrt_x = cbrt_primes.len() as R;
                let pi_sqrt_x = pi_cbrt_x + intermediate_primes.len() as R;
                let constant_term = ((pi_sqrt_x + pi_cbrt_x - 2.0) * (pi_sqrt_x - pi_cbrt_x + 1.0)) / 2.0;
                let mut sqrt_primes_vec = cbrt_primes.to_vec().clone();
                sqrt_primes_vec.extend(intermediate_primes.to_vec().clone());
                let sqrt_primes = Primes::new(sqrt_primes_vec);
                let intermediate_primes_pi_sum = intermediate_primes
                    .to_vec()
                    .iter()
                    .map(|&p| {
                        self.pi_prime(x / p as R, &sqrt_primes)
                    })
                    .sum::<ZPlus>();
                return (legendre_sum + constant_term as Z - intermediate_primes_pi_sum as Z) as ZPlus;
            } else {
                panic!("Not enough in-memory primes to compute pi({})", x);
            }
        };
    }
}

#[cfg(test)]
mod meissel_tests {
    use std::path::Path;

    use crate::prime_table::PrimeTableReader;

    use super::*;

    #[test]
    fn test_meissel_fast() {
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67]);
        let strategy: Meissel = Meissel::new(&primes);
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
    fn test_meissel_slow() {
        if let Some(prime_table_reader) = PrimeTableReader::first_million_from_file() {
            let primes = Primes::new(prime_table_reader.first_million_primes());
            let strategy: Meissel = Meissel::new(&primes);
            assert_eq!(strategy.pi(100000.0_f64), 9592);
        }
    }
}