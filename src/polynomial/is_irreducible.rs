use crate::divisors::divisors::Divisors;
use crate::polynomial::polynomial::{Polynomial, Q, Z};
use crate::prime_factors::prime_factorize::PrimeFactorize;
use crate::primes::primes::Primes;
use num::integer::lcm;
use num::{One, Zero};

impl Polynomial {
    pub fn is_irreducible_over_q(
        &self,
        divisors_strategy: &impl Divisors,
        prime_factorizer: &impl PrimeFactorize,
        primes: &Primes,
    ) -> Option<bool> {
        /*
        If the degree is one, then it is irreducible, because it cannot factor
        into polynomials of lower degree.
        */
        if self.degree() <= 1 {
            return Some(true);
        }

        /*
        If it only has one term, and that term has degree > 0, then it's
        reducible, because 0 will be a root. We just checked the degree, so we
        don't have to check again.
        */
        if self.coefficients.len() == 1 {
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
        From this point onwards, we rely on having integer coefficients, so
        we'll transform the polynomial into one with integer coefficients
        before going on.
        */
        let integer_polynomial = self.convert_to_integer_coefficients();

        /*
        Rational Roots Theorem

        All rational roots of p will have a numerator that divides the constant
        term, and a denominator that divides the leading term coefficient. If
        we take all the combinations of the divisors of the leading and constant
        term coefficients and combine them into a fraction, then pass them
        through the polynomial, then if any value is zero, then clearly the
        polynomial is reducible over the rationals.
        */
        let mut numerators = divisors_strategy
            .divisors(&integer_polynomial.constant_term_coefficient().to_integer());
        numerators.push(Z::from(-1));
        let denominators =
            divisors_strategy.divisors(&integer_polynomial.leading_term_coefficient().to_integer());
        if numerators.iter().any(|n| {
            denominators.iter().any(|d| {
                return self.substitute(Q::new(n.clone(), d.clone())).is_zero();
            })
        }) {
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
        See if we can determine irreducibility using Eisenstein's criterion.
        */
        if Polynomial::is_irreducible_by_eisenstein_criterion(
            &integer_polynomial.try_integer_coefficients(),
            prime_factorizer,
        ) {
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
        */
        let leading_coefficient_prime_factors = prime_factorizer
            .prime_factors(&integer_polynomial.leading_term_coefficient().to_integer());
        let mut primes_vec = primes.to_vec().clone();
        primes_vec.retain(|&q| {
            let q_z_ref = &Z::from(q);
            !leading_coefficient_prime_factors.contains(q_z_ref)
                && !(self.constant_term_coefficient() % q_z_ref).is_zero()
        });
        /*
        Todo: Do something more coherent than arbitrarily taking the first
              five matching primes.
        */
        primes_vec.truncate(5);
        for q in primes_vec {
            let mut coefficients_mod_q = integer_polynomial.try_integer_coefficients();
            for c in coefficients_mod_q.iter_mut() {
                let c_mod_q = c.clone() % Z::from(q);
                *c = c_mod_q;
            }
            if Polynomial::is_irreducible_by_eisenstein_criterion(
                &coefficients_mod_q,
                prime_factorizer,
            ) {
                return Some(true);
            }
        }

        /* Give up and return no answer */
        return None;
    }

    fn is_irreducible_by_eisenstein_criterion(
        coefficients: &Vec<Z>,
        prime_factorizer: &impl PrimeFactorize,
    ) -> bool {
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
        let constant_term = &coefficients[0].clone();
        let leading_coefficient = &coefficients[coefficients.len() - 1].clone();
        let mut non_leading_non_zero_coefficients = coefficients.clone();
        non_leading_non_zero_coefficients.remove(coefficients.len() - 1);
        non_leading_non_zero_coefficients.retain(|f| !f.is_zero());
        let mut common_prime_factors =
            prime_factorizer.prime_factors(&non_leading_non_zero_coefficients[0]);
        for c in non_leading_non_zero_coefficients {
            let prime_factors = prime_factorizer.prime_factors(&c);
            common_prime_factors.retain(|p| prime_factors.contains(p))
        }
        if !common_prime_factors.is_empty() {
            for q in common_prime_factors {
                let q_ref = &q;
                if !(leading_coefficient % q_ref).is_zero()
                    && !(constant_term % &(q_ref * q_ref)).is_zero()
                {
                    return true;
                }
            }
        }
        return false;
    }

    fn convert_to_integer_coefficients(&self) -> Polynomial {
        let mut lcm_accumulator = Z::one();
        for c in &self.coefficients {
            lcm_accumulator = lcm(lcm_accumulator, c.denom().clone());
        }
        let lcm = lcm_accumulator;
        let integer_coefficients: Vec<Z> = self
            .coefficients
            .iter()
            .map(|c| (c * &lcm).to_integer())
            .collect();
        Polynomial::from(integer_coefficients)
    }

    fn try_integer_coefficients(&self) -> Vec<Z> {
        self.coefficients.iter().map(|c| c.to_integer()).collect()
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;
    use crate::divisors::library_divisors::LibraryDivisors;
    use crate::prime_factors::recursive_prime_factorize::RecursivePrimeFactorize;

    fn check_irreducibility(p: Polynomial) -> Option<bool> {
        p.is_irreducible_over_q(
            &LibraryDivisors::new(),
            &RecursivePrimeFactorize::default(),
            &Primes::new(vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]),
        )
    }

    #[test]
    fn test_irreducibility_check_edge_cases() {
        // t + 1
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![1, 1])),
            Some(true)
        );

        // t^5
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![0, 0, 0, 0, 0, 1])),
            Some(false)
        );

        // 10
        assert_eq!(check_irreducibility(Polynomial::from(vec![10])), Some(true));
    }

    #[test]
    fn test_irreducibility_check_rational_roots_theorem() {
        // (t - 3)(t - 2)
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![6, -5, 1])),
            Some(false)
        );

        // (t + 1/2)(t^2 - 5)
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![-5, -10, 1, 2])),
            Some(false)
        );

        // 2t^2 + t + 1
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![1, 1, 2])),
            Some(true)
        );

        // t^3 + 2t^2 - 4
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![1, 1, 2])),
            Some(true)
        );

        // 4t^3 + t^2 - t + 3
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![3, -1, 1, 4])),
            Some(true)
        );

        // 3t^3 + 4t^2 - 6t + 18
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![18, -6, 4, 3])),
            Some(true)
        );
    }

    #[test]
    fn test_irreducibility_check_eisenstein_criterion() {
        // t^4 - 3t + 6
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![6, -3, 0, 0, 1])),
            Some(true)
        );
    }

    #[test]
    fn test_irreducibility_check_taking_mod_p() {
        // 2t^4 + 3t^2 + 3t + 18
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![18, 3, 3, 0, 2])),
            Some(true)
        );
    }

    #[test]
    fn test_irreducibility_no_result() {
        // t^4 + 5t^2 + 4
        assert_eq!(
            check_irreducibility(Polynomial::from(vec![4, 0, 5, 0, 1])),
            None
        );
    }
}
