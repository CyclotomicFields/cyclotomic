use std::cmp::PartialOrd;
use std::fmt::Display;

use crate::prime_counter::prime_counter::PrimeCounter;

type R = f64;
type ZPlus = u64;
type Z = i64;

#[derive(Debug, Eq, PartialEq)]
pub struct Primes {
    primes: Vec<ZPlus>,
}

/*
A collection type for in-memory prime numbers, used to make computations involving bigger prime
numbers. For example, you require a prime number table for fast calculation of large values of
pi(x).
*/
impl Primes {
    pub fn new(primes: Vec<ZPlus>) -> Primes {
        Primes::assert_sorted_ascending(&primes);
        Primes { primes }
    }

    fn new_subset(primes: Vec<ZPlus>) -> Primes {
        Primes { primes }
    }

    fn assert_sorted_ascending<T: PartialOrd + Display>(vec: &Vec<T>) {
        assert!(!vec.is_empty());
        let mut current_max: &T = &vec[0];
        for j in 0..vec.len() {
            assert!(
                vec[j] >= *current_max,
                "{} was not greater than or equal to a previous \
                    element in the list of primes, {}",
                vec[j],
                *current_max
            );
            current_max = &vec[j];
        }
    }

    pub fn to_vec(&self) -> &Vec<ZPlus> {
        return &self.primes;
    }

    /*
    Returns the first n primes in this list
    */
    pub fn first_n(&self, n: ZPlus) -> Option<Primes> {
        if n > self.primes.len() as ZPlus {
            eprintln!("Don't have {} primes", n);
            return None;
        }
        if n == 0 {
            eprintln!("Unexpectedly asked for 0 primes");
            return None;
        }
        let mut primes_clone = self.primes.clone();
        primes_clone.truncate(n as usize);
        return Some(Primes::new_subset(primes_clone));
    }

    /*
    Returns the number of primes not greater than x.
    */
    fn pi_prime(&self, x: ZPlus) -> Option<ZPlus> {
        if x > self.primes[self.primes.len() - 1] {
            eprintln!(
                "Can't evaluate pi({}) with primes in the interval ({}, {})",
                x,
                self.primes[0],
                self.primes[self.primes.len() - 1]
            );
            return None;
        } else if x == self.primes[self.primes.len() - 1] {
            return Some(self.primes.len() as ZPlus);
        }

        let mut i = 0;
        let mut p = self.primes[i];
        loop {
            if p > x {
                return Some(i as ZPlus);
            }
            if i == self.primes.len() - 1 {
                eprintln!("Ran out of primes while evaluating pi({})", x);
                return None;
            }
            i += 1;
            p = self.primes[i];
        }
    }

    /*
    Returns all the primes greater than the lower bound, and not greater than the upper bound.
    */
    pub fn range(&self, lower_bound: ZPlus, upper_bound: ZPlus) -> Option<Primes> {
        if upper_bound > self.primes[self.primes.len() - 1] || upper_bound < lower_bound {
            eprintln!(
                "Can't get primes in the range {} to {} using the given {} primes",
                lower_bound,
                upper_bound,
                self.primes.len()
            );
            return None;
        } else if upper_bound == lower_bound {
            return Some(Primes::new_subset(vec![]));
        }
        let mut primes_clone = self.primes.clone();
        primes_clone.retain(|&p| lower_bound < p && p <= upper_bound);
        return Some(Primes::new_subset(primes_clone));
    }

    /*
    Returns the number of primes greater than the lower bound, and not greater than the upper bound.
    */
    pub fn pi_range(&self, lower_bound: ZPlus, upper_bound: ZPlus) -> Option<ZPlus> {
        if upper_bound < lower_bound {
            eprintln!(
                "Upper bound {} is less than lower bound {}",
                upper_bound, lower_bound
            );
            return None;
        }
        let lower_pi: Option<u64> = self.pi_prime(lower_bound);
        if lower_pi.is_none() {
            eprintln!("Couldn't get pi({})", lower_bound);
            return None;
        }
        let upper_pi: Option<u64> = self.pi_prime(upper_bound);
        if upper_pi.is_none() {
            eprintln!("Couldn't get pi({})", upper_bound);
            return None;
        }
        return Some(upper_pi.unwrap() - lower_pi.unwrap());
    }

    pub fn len(&self) -> usize {
        return self.primes.len();
    }

    pub fn nth(&self, x: usize) -> Option<ZPlus> {
        if x == 0 {
            return None;
        }
        return Some(self.primes[x - 1]);
    }
}

impl PrimeCounter for Primes {
    fn pi(&self, x: f64) -> u64 {
        if let Some(pi) = self.pi_prime(x as ZPlus) {
            return pi;
        } else {
            panic!("Didn't have enough primes to evaluate pi({})", x);
        }
    }
}

#[cfg(test)]
pub mod primes_test {
    use std::path::Path;

    use super::*;

    use crate::primes::prime_table::PrimeTableReader;

    #[test]
    fn test_first_n() {
        let primes_vec = vec![2, 3, 5, 7, 11, 13];
        let primes = Primes::new(primes_vec);
        assert_eq!(primes.first_n(0), None);
        assert_eq!(*primes.first_n(1).unwrap().to_vec(), vec![2]);
        assert_eq!(*primes.first_n(2).unwrap().to_vec(), vec![2, 3]);
        assert_eq!(*primes.first_n(3).unwrap().to_vec(), vec![2, 3, 5]);
        assert_eq!(*primes.first_n(4).unwrap().to_vec(), vec![2, 3, 5, 7]);
        assert_eq!(*primes.first_n(5).unwrap().to_vec(), vec![2, 3, 5, 7, 11]);
        assert_eq!(
            *primes.first_n(6).unwrap().to_vec(),
            vec![2, 3, 5, 7, 11, 13]
        );
    }

    #[test]
    fn test_pi() {
        let primes_vec = vec![2, 3, 5, 7, 11, 13];
        let primes = Primes::new(primes_vec);
        assert_eq!(primes.pi_prime(1), Some(0));
        assert_eq!(primes.pi_prime(2), Some(1));
        assert_eq!(primes.pi_prime(3), Some(2));
        assert_eq!(primes.pi_prime(10), Some(4));
        assert_eq!(primes.pi_prime(13), Some(6));
        assert_eq!(primes.pi_prime(15), None);
    }

    #[test]
    fn test_range() {
        let primes_vec = vec![2, 3, 5, 7, 11, 13];
        let primes = Primes::new(primes_vec);
        assert_eq!(*primes.range(0, 5).unwrap().to_vec(), vec![2, 3, 5]);
        assert_eq!(*primes.range(3, 5).unwrap().to_vec(), vec![5]);
        assert_eq!(*primes.range(7, 13).unwrap().to_vec(), vec![11, 13]);
        assert_eq!(*primes.range(2, 2).unwrap().to_vec(), vec![] as Vec<u64>);
        assert_eq!(primes.range(0, 15), None);
        assert_eq!(primes.range(15, 20), None);
        assert_eq!(primes.range(10, 20), None);
        assert_eq!(primes.range(13, 2), None);
    }

    #[test]
    fn test_pi_range() {
        let primes_vec = vec![2, 3, 5, 7, 11, 13];
        let primes = Primes::new(primes_vec);
        assert_eq!(primes.pi_range(0, 5), Some(3));
        assert_eq!(primes.pi_range(3, 5), Some(1));
        assert_eq!(primes.pi_range(7, 13), Some(2));
        assert_eq!(primes.pi_range(5, 5), Some(0));
        assert_eq!(primes.pi_range(15, 20), None);
        assert_eq!(primes.pi_range(10, 20), None);
        assert_eq!(primes.pi_range(13, 2), None);
    }

    #[test]
    fn test_primes_pi_slow() {
        if let Some(prime_table_reader) = PrimeTableReader::first_million_from_file() {
            let strategy = Primes::new(prime_table_reader.first_million_primes());
            assert_eq!(strategy.pi(1000000.0_f64), 78498);
        }
    }
}
