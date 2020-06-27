extern crate num;

use self::num::{BigInt, BigRational, Zero, ToPrimitive, One, Integer, Signed};
use std::ops::{Div, Rem, Sub};
use crate::polynomial::polynomial::Polynomial;
use std::cmp::max;
use std::borrow::Borrow;

type Z = num::bigint::BigInt;

pub struct Euclid {}

impl Euclid {
    pub fn new() -> Euclid {
        Euclid {}
    }

    fn extended_euclidean_augmented_matrix_result_rows(&self, m: &Z, n: &Z) -> (Vec<Z>, Vec<Z>) {
        /*
        Extended Euclidean Algorithm (EEA)

        Bézout's identity states that the GCD of any two numbers can be
        represented as a linear combination of those numbers. We exploit use
        this property in the algorithm to represent our progress in Euclid's
        algorithm towards determining the GCD.
        That is, m*k1 + n*k2 = gcd(m, n) for some integers k1 and k2.

        The algorithm proceeds by constructing a matrix and augmenting it by
        logically appending a row to the bottom, in each iteration. The initial
        matrix represents the a trivial pair of linear combinations for m and n
        and it looks like this, where m >= n.

        m 1 0
        n 0 1

        Every row represents a linear combination of m and n. Each iteration of
        the algorithm will result in appending a new row to the bottom of this
        matrix. To construct the next row, we would take the quotient and of a
        division m / n. The row is then the second-to-bottom row of the matrix
        minus the quotient multiplied by the bottom row of the matrix. In this
        case, that would be (m 1 0) - k * (n 0 1) where k is the integer part
        of m / n.

        We repeat this process, continually adding rows to the bottom of the
        matrix, until we get a remainder of 0, that is, the value in the
        left-most column of the bottom row is 0. At this point, we can say that
        the GCD of n and m is the remainder in the second-to-bottom row,
        according to Euclid's original algorithm.
        */
        let mut trailing_row = vec![(*m).clone(), Z::one(), Z::zero()];
        let mut leading_row = vec![(*n).clone(), Z::zero(), Z::one()];
        while !&leading_row[0].is_zero() {
            let (quotient, remainder) = (&trailing_row[0]).div_rem(&leading_row[0]);
            let new_leading_row = vec![remainder,
                                       &trailing_row[1] - &quotient * &leading_row[1],
                                       &trailing_row[2] - &quotient * &leading_row[2]];
            trailing_row = leading_row;
            leading_row = new_leading_row;
        }
        /*
        These two vectors represent a pair of identities of the form:

        gcd(m, n) = m * k1 + n * k2
        0         = m * k3 + n * k4

        Where we are returning just the coefficients:

        (gcd(m, n) k1 k2)
        (0         k3 k4)
        */
        (trailing_row, leading_row)
    }

    pub fn gcd(&self, m: &Z, n: &Z) -> Z {
        if m == n {
            return (*m).clone();
        } else if n > m {
            return self.gcd(n, m);
        }
        let (gcd_identity, _) = self.extended_euclidean_augmented_matrix_result_rows(m, n);
        return gcd_identity[0].clone();
    }

    pub fn multiplicative_inverse_mod_n(&self, i: &Z, n: &Z) -> Option<Z> {
        let (gcd_identity, _) = self.extended_euclidean_augmented_matrix_result_rows(n, i);
        if !gcd_identity[0].is_one() {
            /*
            If the GCD isn't one, then it means that i and n aren't relatively
            prime, which means that i has no multiplicative in the
            integers mod n. In this case, we return None to indicate that there
            is no multiplicative inverse for this value in this ring.
            */
            None
        } else {
            /*
            If there is a multiplicative inverse, then the EEA iteration has
            given us an identity 1 = n * k1 + i * k2 for integers k1 and k2.
            This is represented in coefficient form as (1 k1 k2). Since we are
            in the integers mod n, the n * k1 component is ignored. To get the
            multiplicative inverse of i, we take the coefficient k2 modulo n.
            */
            Some(gcd_identity[2].mod_floor(n))
        }
    }

    pub fn multiplicative_inverse_mod_polynomial(&self, p: Polynomial, m: Polynomial) -> Option<Polynomial> {
        /*
        As with finding the multiplicative inverse for integers mod n, we
        take the GCD of the polynomial with the principal ideal polynomial of
        the polynomial ring, and then find an associated value within the ring.

        I've copied the implementation from above because I can't be bothered to
        do the type-fu required to make it more generic.
        */
        let mut trailing_row = vec![m, Polynomial::one(), Polynomial::zero()];
        let mut leading_row = vec![p, Polynomial::zero(), Polynomial::one()];
        while !(leading_row[0].is_zero()) {
            let (quotient, remainder) = (&trailing_row[0]).div(&leading_row[0]);
            let new_leading_row = vec![remainder,
                                       (&trailing_row[1]).sub(&leading_row[1].mul(&quotient)),
                                       (&trailing_row[2]).sub(&leading_row[2].mul(&quotient))];
            trailing_row = leading_row;
            leading_row = new_leading_row;
        }
        if trailing_row[0].is_one() {
            Some(trailing_row[2].clone())
        } else if trailing_row[0].neg().is_one() {
            Some(trailing_row[2].clone().neg())
        } else {
            None
        }
    }
}

#[cfg(test)]
mod euclid_tests {
    use super::*;

    #[test]
    fn test_multiplicative_inverse_mod_polynomial() {
        let euclid = Euclid::new();
        // (t + 2)^-1 == t^2 + t + 1 in the ring Z[t] / t^3 + 3t^2 + 3t + 1
        assert_eq!(euclid.multiplicative_inverse_mod_polynomial(
            Polynomial::from(vec![2, 1]), Polynomial::from(vec![1, 3, 3, 1])),
                   Some(Polynomial::from(vec![1, 1, 1])));

        // (t^3 - 1)^-1 has no inverse in the ring Z[t] / t^5 + 4t^3 + 2
        assert_eq!(euclid.multiplicative_inverse_mod_polynomial(
            Polynomial::from(vec![-1, 0, 0, 1]), Polynomial::from(vec![2, 0, 0, 4, 0, 1])),
                   None);

        // (2t  + 4)^-1 == (-1/2)t^3 - t^2 - t in the ring Z[t] / t^4 + 4t^3 + 6t^2 + 4t + 1
        assert_eq!(euclid.multiplicative_inverse_mod_polynomial(
            Polynomial::from(vec![4, 2]), Polynomial::from(vec![1, 4, 6, 4, 1])),
                   Some(Polynomial::from_small_fractions(vec![0, -1, -1, -1], vec![1, 1, 1, 2])));
    }

    #[test]
    fn test_extended_euclid_gcd() {
        let euclid = Euclid::new();
        assert_eq!(euclid.gcd(&Z::from(80), &Z::from(62)), Z::from(2));
        assert_eq!(euclid.gcd(&Z::from(5), &Z::from(5)), Z::from(5));
        assert_eq!(euclid.gcd(&Z::from(420), &Z::from(1782)), Z::from(6));
        assert_eq!(euclid.gcd(&Z::from(371250), &Z::from(873630)), Z::from(90));
        assert_eq!(euclid.gcd(&Z::from(3283741336253_i64), &Z::from(0234551825263_i64)), Z::from(7));
    }

    #[test]
    fn test_multiplicative_inverse_mod_n() {
        let euclid = Euclid::new();
        assert_eq!(euclid.multiplicative_inverse_mod_n(&Z::from(5), &Z::from(17)), Some(Z::from(7)));
        assert_eq!(euclid.multiplicative_inverse_mod_n(&Z::from(3), &Z::from(26)), Some(Z::from(9)));
        assert_eq!(euclid.multiplicative_inverse_mod_n(&Z::from(5), &Z::from(10)), None);
        assert_eq!(euclid.multiplicative_inverse_mod_n(&Z::from(16), &Z::from(24)), None);
        assert_eq!(euclid.multiplicative_inverse_mod_n(&Z::from(27), &Z::from(392)), Some(Z::from(363)));
    }
}