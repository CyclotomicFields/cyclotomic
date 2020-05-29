extern crate num;

use num::pow::pow;
use std::cmp::PartialOrd;
use std::fmt::Display;
use self::num::BigInt;

type Z = num::bigint::BigInt;
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
        let mut sum: Z = BigInt::from(0);
        for j in 0..self.coefficients.len() {
            sum += self.coefficients[j].clone() * num::pow(t.clone(), self.degrees[j]);
        }
        return sum;
    }

    pub fn is_monic(&self) -> bool {
        return self.coefficients[self.coefficients.len() - 1] == BigInt::from(1);
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
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    fn vec_z(vec: Vec<i64>) -> Vec<Z> {
        vec.iter().map(|&i| BigInt::from(i)).collect()
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

        assert_eq!(p.substitute(BigInt::from(5)), BigInt::from(1));
        assert_eq!(p.substitute(BigInt::from(0)), BigInt::from(1));
        assert_eq!(p.substitute(BigInt::from(-1)), BigInt::from(1));
        assert_eq!(p.substitute(BigInt::from(-14)), BigInt::from(1));
        assert_eq!(p.substitute(BigInt::from(1)), BigInt::from(1));
        assert_eq!(p.substitute(BigInt::from(2)), BigInt::from(1));
        assert_eq!(p.substitute(BigInt::from(125716)), BigInt::from(1));

        // p = 2t + 1
        p = Polynomial::new(vec_z(vec![1, 2]), vec![0, 1]);

        assert_eq!(p.substitute(BigInt::from(5)), BigInt::from(11));
        assert_eq!(p.substitute(BigInt::from(0)), BigInt::from(1));
        assert_eq!(p.substitute(BigInt::from(-1)), BigInt::from(-1));
        assert_eq!(p.substitute(BigInt::from(-14)), BigInt::from(-27));
        assert_eq!(p.substitute(BigInt::from(1)), BigInt::from(3));
        assert_eq!(p.substitute(BigInt::from(2)), BigInt::from(5));
        assert_eq!(p.substitute(BigInt::from(125716)), BigInt::from(251433));

        // p = t^5 - t^3 - 2t^2 - 1
        p = Polynomial::new(vec_z(vec![-1, 0, -2, -1, 0, 1]), vec![0, 1, 2, 3, 4, 5]);

        assert_eq!(p.substitute(BigInt::from(5)), BigInt::from(2949));
        assert_eq!(p.substitute(BigInt::from(0)), BigInt::from(-1));
        assert_eq!(p.substitute(BigInt::from(-1)), BigInt::from(-3));
        assert_eq!(p.substitute(BigInt::from(-14)), BigInt::from(-535473));
        assert_eq!(p.substitute(BigInt::from(1)), BigInt::from(-3));
        assert_eq!(p.substitute(BigInt::from(2)), BigInt::from(15));
        assert_eq!(p.substitute(BigInt::from(123)), BigInt::from(28151165717_i64));
    }
}