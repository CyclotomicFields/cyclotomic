extern crate num;

use num::pow::pow;
use std::cmp::PartialOrd;
use std::fmt::Display;

type Z = i64;
type ZPlus = usize;

pub struct Polynomial {
    coefficients: Vec<Z>,
    degrees: Vec<ZPlus>,
}

impl Polynomial {
    pub fn new(coefficients: Vec<Z>, degrees: Vec<ZPlus>) -> Polynomial {
        assert_eq!(coefficients.len(), degrees.len());
        Polynomial::assert_sorted_ascending(&degrees);
        return Polynomial { coefficients, degrees };
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

    pub fn substitute(&self, t: Z) -> Z {
        let mut sum: Z = 0;
        for j in 0..self.coefficients.len() {
            sum += self.coefficients[j] * pow(t, self.degrees[j]);
        }
        return sum;
    }

    pub fn is_monic(&self) -> bool {
        return self.coefficients[self.coefficients.len() - 1] == 1;
    }

    pub fn is_irreducible_over_z() -> Option<bool> {
        // First check: Rational Roots Theorem
        // All rational roots of p will have a numerator that divides the constant term, and a
        // denominator that divides the leading term coefficient. If we take all the combinations of
        // the divisors of the leading and constant term coefficients and combine them into a
        // fraction, then pass them through the polynomial, then if any value is zero, then clearly
        // the polynomial is reducible over the rationals.

        // Second check: Eisenstein's Criterion
        // If there exists a prime, q, which is a factor of every non-leading term, not a factor of
        // the leading term, and also where q squared is not a factor of the constant term, then
        // the polynomial is irreducible.

        // At this point we give up and say "I don't know".
        return None;
    }

    fn divisors(i: Z) -> Vec<Z> {
        return vec![];
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_monic_check() {
        // p = 1
        assert!(Polynomial::new(vec![1], vec![0]).is_monic());
        // p = 2t + 1
        assert!(!Polynomial::new(vec![1, 2], vec![0, 1]).is_monic());
        // p = t^5 - t^3 - 2t^2 - 1
        assert!(Polynomial::new(vec![-1, 0, -2, -1, 0, 1], vec![0, 1, 2, 3, 4, 5]).is_monic());
    }

    #[test]
    fn test_polynomial_substitution() {
        // p = 1
        let mut p: Polynomial = Polynomial::new(vec![1], vec![0]);

        assert_eq!(p.substitute(5), 1);
        assert_eq!(p.substitute(0), 1);
        assert_eq!(p.substitute(-1), 1);
        assert_eq!(p.substitute(-14), 1);
        assert_eq!(p.substitute(1), 1);
        assert_eq!(p.substitute(2), 1);
        assert_eq!(p.substitute(125716), 1);

        // p = 2t + 1
        p = Polynomial::new(vec![1, 2], vec![0, 1]);

        assert_eq!(p.substitute(5), 11);
        assert_eq!(p.substitute(0), 1);
        assert_eq!(p.substitute(-1), -1);
        assert_eq!(p.substitute(-14), -27);
        assert_eq!(p.substitute(1), 3);
        assert_eq!(p.substitute(2), 5);
        assert_eq!(p.substitute(125716), 251433);

        // p = t^5 - t^3 - 2t^2 - 1
        p = Polynomial::new(vec![-1, 0, -2, -1, 0, 1], vec![0, 1, 2, 3, 4, 5]);

        assert_eq!(p.substitute(5), 2949);
        assert_eq!(p.substitute(0), -1);
        assert_eq!(p.substitute(-1), -3);
        assert_eq!(p.substitute(-14), -535473);
        assert_eq!(p.substitute(1), -3);
        assert_eq!(p.substitute(2), 15);
        assert_eq!(p.substitute(123), 28151165717);
    }
}