extern crate num;

use self::num::{BigInt, BigRational, Zero, ToPrimitive, One, Integer};
use std::ops::{Div, Rem};

type Z = num::bigint::BigInt;

pub struct Euclid {}

impl Euclid {
    pub fn new() -> Euclid {
        Euclid {}
    }

    pub fn gcd(&self, m: &Z, n: &Z) -> Z {
        if m == n {
            return (*m).clone();
        } else if n > m {
            return self.gcd(n, m);
        }
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
        return trailing_row[0].clone();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extended_euclid_gcd() {
        let euclid = Euclid::new();
        assert_eq!(euclid.gcd(&Z::from(80), &Z::from(62)), Z::from(2));
        assert_eq!(euclid.gcd(&Z::from(5), &Z::from(5)), Z::from(5));
        assert_eq!(euclid.gcd(&Z::from(420), &Z::from(1782)), Z::from(6));
        assert_eq!(euclid.gcd(&Z::from(371250), &Z::from(873630)), Z::from(90));
        assert_eq!(euclid.gcd(&Z::from(3283741336253_i64), &Z::from(0234551825263_i64)), Z::from(7));
    }
}