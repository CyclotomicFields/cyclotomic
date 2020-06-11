extern crate num;

use num::pow::pow;
use std::cmp::PartialOrd;
use std::fmt::Display;
use self::num::{BigInt, BigRational, Zero, ToPrimitive, Integer, One};
use crate::divisors::divisors;
use crate::divisors::library_divisors::LibraryDivisors;
use crate::divisors::divisors::Divisors;
use crate::prime_factors::recursive_prime_factorize::RecursivePrimeFactorize;
use crate::prime_factors::prime_factorize::PrimeFactorize;
use crate::primes::primes::Primes;
use std::ops::Mul;

type Z = num::bigint::BigInt;
type Q = num::rational::BigRational;
type ZPlus = usize;

pub struct Polynomial {
    coefficients: Vec<Z>,
    degrees: Vec<ZPlus>,
}

impl Polynomial {
    pub fn new(coefficients: Vec<Z>, degrees: Vec<ZPlus>) -> Polynomial {
        assert_eq!(coefficients.len(), degrees.len());
        Polynomial::assert_sorted_ascending(&degrees);
        return Polynomial {
            coefficients,
            degrees,
        };
    }

    fn assert_sorted_ascending<T: PartialOrd + Display>(vec: &Vec<T>) {
        let mut current_max: &T = &vec[0];
        for j in 0..vec.len() {
            assert!(vec[j] >= *current_max,
                    "{} was not greater than or equal to a previous \
                    element in the list of degrees, {}", vec[j], *current_max);
            current_max = &vec[j];
        }
    }

    pub fn substitute(&self, t: Q) -> Q {
        let mut sum: Q = Q::zero();
        for j in 0..self.coefficients.len() {
            sum += num::pow(t.clone(), self.degrees[j]).mul(&self.coefficients[j]);
        }
        return sum;
    }

    pub fn is_monic(&self) -> bool {
        return self.leading_term_coefficient().is_one();
    }

    pub fn degree(&self) -> ZPlus {
        return self.degrees[self.degrees.len() - 1];
    }

    pub fn leading_term_coefficient(&self) -> &Z {
        &self.coefficients[self.coefficients.len() - 1]
    }

    pub fn constant_term_coefficient(&self) -> &Z {
        &self.coefficients[0]
    }

    pub fn is_irreducible_over_q(&self) -> Option<bool> {
        /*
        If the degree is one, then it is irreducible, because it cannot factor
        into polynomials of lower degree.
        */
        if self.degree() <= 1 {
            return Some(true);
        }

        /*
        If it only has one term, and that term has degree >= 2, then it's
        reducible, because 0 will be a root.
        */
        if self.degrees.len() == 1 && self.degree() >= 2 {
            return Some(false);
        }

        /*
        If the constant term is 0, then it's reducible, because 0 will be a
        root.
        */
        if self.constant_term_coefficient().is_zero() {
            return Some(false);
        }

        /*
        Rational Roots Theorem

        All rational roots of p will have a numerator that divides the constant
        term, and a denominator that divides the leading term coefficient. If
        we take all the combinations of the divisors of the leading and constant
        term coefficients and combine them into a fraction, then pass them
        through the polynomial, then if any value is zero, then clearly the
        polynomial is reducible over the rationals.
        */
        let divisors_strategy = LibraryDivisors::new();
        let mut numerators = divisors_strategy.divisors(self.constant_term_coefficient());
        numerators.push(Z::from(-1));
        let denominators = divisors_strategy.divisors(self.leading_term_coefficient());
        if numerators.iter().any(|n| denominators.iter().any(|d| {
            return self.substitute(Q::new(n.clone(), d.clone())).is_zero();
        })) {
            return Some(false);
        }

        /*
        All factorisations of degree 2 or degree 3 polynomials must result in a
        degree 1 factor, also known as a root. Therefore if a polynomial of
        degree 2 or degree 3 has no roots, then there are no factorisations,
        which means that it is an irreducible polynomial.π
        */
        if self.degree() <= 3 {
            return Some(true);
        }

        /*
        Eisenstein's Criterion

        The polynomial is irreducible if there exists a prime, q, such that:
          - q is a factor of every non-leading term
          - q is not a factor of the leading term
          - q squared is not a factor of the constant term

        Find the common prime factors of all the non-leading coefficients, and
        check all of them against the second and third points of the criterion,
        as listed above.
        */
        let prime_factorizer = RecursivePrimeFactorize::default();
        if self.is_irreducible_by_eisenstein_criterion(&self.coefficients, &prime_factorizer) {
            return Some(true);
        }

        /*
        Reducing mod q

        If q is a prime that is not a factor of the leading coefficient, then
        if the polynomial is irreducible over Z mod q, then it is also
        irreducible over Q.

        We can check this for some sensible number of candidate primes, q. When
        checking for irreducibility over Z mod q, we use the Eisenstein
        criterion. We aren't interested in the case where the polynomial is
        reducible, because that isn't conclusive information.
        Todo: It may be computationally faster to use the rational roots theorem
              to fail fast for reducible polynomials.
        Todo: Make the model coherent. Polynomial instances shouldn't
              have their own instance of Primes, they should share a big
              blob of central, read-only primes. Same for the other objects
              created in this method.
        */
        let primes = Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
        let leading_coefficient_prime_factors = prime_factorizer.prime_factors(self.leading_term_coefficient());
        let mut primes_vec = primes.to_vec().clone();
        primes_vec.retain(|&q| {
            let q_z_ref = &Z::from(q);
            !leading_coefficient_prime_factors.contains(q_z_ref)
                && !(self.constant_term_coefficient() % q_z_ref).is_zero()
        });
        /*
        Todo: Do something more coherent than arbitrarily taking the first
              five matching primes out of an already arbitrarily chosen list
              of primes.
        */
        primes_vec.truncate(5);
        for q in primes_vec {
            let mut coefficients_mod_q = self.coefficients.clone();
            for c in coefficients_mod_q.iter_mut() {
                let c_mod_q = c.clone() % Z::from(q);
                *c = c_mod_q;
            }
            if self.is_irreducible_by_eisenstein_criterion(&coefficients_mod_q, &prime_factorizer) {
                return Some(true);
            }
        }

        /* Give up and return no answer */
        return None;
    }

    fn is_irreducible_by_eisenstein_criterion(&self, coefficients: &Vec<Z>, prime_factorizer: &impl PrimeFactorize) -> bool {
        let constant_term = &coefficients[0].clone();
        let leading_coefficient = &coefficients[coefficients.len() - 1].clone();
        let mut non_leading_non_zero_coefficients = coefficients.clone();
        non_leading_non_zero_coefficients.remove(coefficients.len() - 1);
        non_leading_non_zero_coefficients.retain(|f| !f.is_zero());
        let mut common_prime_factors = prime_factorizer.prime_factors(&non_leading_non_zero_coefficients[0]);
        for c in non_leading_non_zero_coefficients {
            let prime_factors = prime_factorizer.prime_factors(&c);
            common_prime_factors.retain(|p| prime_factors.contains(p))
        }
        if !common_prime_factors.is_empty() {
            for q in common_prime_factors {
                let q_ref = &q;
                if !(leading_coefficient % q_ref).is_zero()
                    && !(constant_term % &(q_ref * q_ref)).is_zero() {
                    return true;
                }
            }
        }
        return false;
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    fn vec_z(vec: Vec<i64>) -> Vec<Z> {
        vec.iter().map(|&i| BigInt::from(i)).collect()
    }

    fn q_from_i64(n: i64) -> Q {
        return Q::from(Z::from(n));
    }

    #[test]
    fn test_irreducibility_check_edge_cases() {
        // t + 1
        assert_eq!(Polynomial::new(vec_z(vec![1, 1]), vec![0, 1]).is_irreducible_over_q(), Some(true));

        // t^5
        assert_eq!(Polynomial::new(vec_z(vec![0, 0, 0, 0, 0, 1]), vec![0, 1, 2, 3, 4, 5]).is_irreducible_over_q(), Some(false));

        // 10
        assert_eq!(Polynomial::new(vec_z(vec![10]), vec![0]).is_irreducible_over_q(), Some(true));
    }

    #[test]
    fn test_irreducibility_check_rational_roots_theorem() {
        // (t - 3)(t - 2)
        assert_eq!(Polynomial::new(vec_z(vec![6, -5, 1]), vec![0, 1, 2]).is_irreducible_over_q(), Some(false));

        // (t + 1/2)(t^2 - 5)
        assert_eq!(Polynomial::new(vec_z(vec![-5, -10, 1, 2]), vec![0, 1, 2, 3]).is_irreducible_over_q(), Some(false));

        // 2t^2 + t + 1
        assert_eq!(Polynomial::new(vec_z(vec![1, 1, 2]), vec![0, 1, 2]).is_irreducible_over_q(), Some(true));

        // t^3 + 2t^2 - 4
        assert_eq!(Polynomial::new(vec_z(vec![1, 1, 2]), vec![0, 1, 2]).is_irreducible_over_q(), Some(true));

        // 4t^3 + t^2 - t + 3
        assert_eq!(Polynomial::new(vec_z(vec![3, -1, 1, 4]), vec![0, 1, 2, 3]).is_irreducible_over_q(), Some(true));

        // 3t^3 + 4t^2 - 6t + 18
        assert_eq!(Polynomial::new(vec_z(vec![18, -6, 4, 3]), vec![0, 1, 2, 3]).is_irreducible_over_q(), Some(true));
    }

    #[test]
    fn test_irreducibility_check_eisenstein_criterion() {
        // t^4 - 3t + 6
        assert_eq!(Polynomial::new(vec_z(vec![6, -3, 0, 0, 1]), vec![0, 1, 2, 3, 4]).is_irreducible_over_q(), Some(true));
    }

    #[test]
    fn test_irreducibility_check_taking_mod_p() {
        // 2t^4 + 3t^2 + 3t + 18
        assert_eq!(Polynomial::new(vec_z(vec![18, 3, 3, 0, 2]), vec![0, 1, 2, 3, 4]).is_irreducible_over_q(), Some(true));
    }

    #[test]
    fn test_irreducibility_no_result() {
        // t^4 + 5t^2 + 4
        assert_eq!(Polynomial::new(vec_z(vec![4, 0, 5, 0, 1]), vec![0, 1, 2, 3, 4]).is_irreducible_over_q(), None);
    }

    #[test]
    fn test_monic_check() {
        // p = 1
        assert!(Polynomial::new(vec_z(vec![1]), vec![0]).is_monic());
        // p = 2t + 1
        assert!(!Polynomial::new(vec_z(vec![1, 2]), vec![0, 1]).is_monic());
        // p = t^5 - t^3 - 2t^2 - 1
        assert!(Polynomial::new(vec_z(vec![-1, 0, -2, -1, 0, 1]), vec![0, 1, 2, 3, 4, 5]).is_monic());
    }

    #[test]
    fn test_polynomial_substitution() {
        // p = 1
        let mut p: Polynomial = Polynomial::new(vec_z(vec![1]), vec![0]);

        assert_eq!(p.substitute(q_from_i64(5)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(0)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(-1)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(-14)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(1)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(2)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(125716)), q_from_i64(1));

        // p = 2t + 1
        p = Polynomial::new(vec_z(vec![1, 2]), vec![0, 1]);

        assert_eq!(p.substitute(q_from_i64(5)), q_from_i64(11));
        assert_eq!(p.substitute(q_from_i64(0)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(-1)), q_from_i64(-1));
        assert_eq!(p.substitute(q_from_i64(-14)), q_from_i64(-27));
        assert_eq!(p.substitute(q_from_i64(1)), q_from_i64(3));
        assert_eq!(p.substitute(q_from_i64(2)), q_from_i64(5));
        assert_eq!(p.substitute(q_from_i64(125716)), q_from_i64(251433));

        // p = t^5 - t^3 - 2t^2 - 1
        p = Polynomial::new(vec_z(vec![-1, 0, -2, -1, 0, 1]), vec![0, 1, 2, 3, 4, 5]);

        assert_eq!(p.substitute(q_from_i64(5)), q_from_i64(2949));
        assert_eq!(p.substitute(q_from_i64(0)), q_from_i64(-1));
        assert_eq!(p.substitute(q_from_i64(-1)), q_from_i64(-3));
        assert_eq!(p.substitute(q_from_i64(-14)), q_from_i64(-535473));
        assert_eq!(p.substitute(q_from_i64(1)), q_from_i64(-3));
        assert_eq!(p.substitute(q_from_i64(2)), q_from_i64(15));
        assert_eq!(p.substitute(q_from_i64(123)), q_from_i64(28151165717_i64));
    }
}