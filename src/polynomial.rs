extern crate num;
extern crate divisors;

use num::pow::pow;
use std::cmp::PartialOrd;
use std::fmt::Display;
use self::num::{BigInt, BigRational, Zero, ToPrimitive};
use self::divisors::get_divisors;

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

    pub fn substitute(&self, t: Q) -> Q {
        let mut sum: Q = Q::zero();
        for j in 0..self.coefficients.len() {
            sum += Q::from(self.coefficients[j].clone()) * num::pow(t.clone(), self.degrees[j]);
        }
        return sum;
    }

    pub fn is_monic(&self) -> bool {
        return self.coefficients[self.coefficients.len() - 1] == BigInt::from(1);
    }

    pub fn is_irreducible_over_z(&self) -> Option<bool> {
        /*
        First check: Rational Roots Theorem

        All rational roots of p will have a numerator that divides the constant
        term, and a denominator that divides the leading term coefficient. If
        we take all the combinations of the divisors of the leading and constant
        term coefficients and combine them into a fraction, then pass them
        through the polynomial, then if any value is zero, then clearly the
        polynomial is reducible over the rationals.
        */
        let numerators = Polynomial::divisors(self.coefficients[0].clone());
        let denominators = Polynomial::divisors(self.coefficients[self.coefficients.len() - 1].clone());
        if numerators.iter().any(|n| denominators.iter().any(|d| {
            return self.substitute(Q::new(n.clone(), d.clone())).is_zero();
        })) {
            return Some(false);
        }

        /*
        Second check: Eisenstein's Criterion

        The polynomial is irreducible if there exists a prime, q, such that:
          - q is a factor of every non-leading term
          - q is not a factor of the leading term
          - q squared is not a factor of the constant term

        */

        /* Give up and return no answer */
        return None;
    }

    /*
    Returns a vector containing all the divisors of z, including 1 and itself.
    */
    fn divisors(z: Z) -> Vec<Z> {
        let z_u128 = z.to_u128().unwrap();
        let mut divisors = get_divisors(z_u128);
        divisors.push(1);
        if z_u128 != 1 {
            divisors.push(z_u128);
        }
        return divisors.iter().map(|&u| Z::from(u)).collect();
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;
    use super::divisors::get_divisors;

    fn vec_z(vec: Vec<i64>) -> Vec<Z> {
        vec.iter().map(|&i| BigInt::from(i)).collect()
    }

    fn q_from_i64(n: i64) -> Q {
        return Q::from(Z::from(n));
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

    #[test]
    fn test_irreducibility_check_rational_roots_theorem() {
        let p = Polynomial::new(vec_z(vec![6, -5, 1]), vec![0, 1, 2]);
        assert_eq!(p.is_irreducible_over_z(), Some(false))
    }
}