use std::ops::{Mul, MulAssign};

use num::{One, one, Zero};

use crate::polynomial::polynomial::{Polynomial, Z, Q};

impl Polynomial {
    pub fn mul_mut(&mut self, rhs: &Self) {
        /*
        Matrix product method

        Represent the coefficients of the lhs polynomial as diagonal entries in
        a not-necessarily-square matrix, and perform a matrix multiplication
        with a column vector containing the coefficients of the rhs polynomial.

        I have no idea if this is fast. I thought it might be, if we exploit a
        fast library for doing matrix multiplications, which is smart enough to
        exploit vector instructions on the CPU. Even then I'd have to code it
        in a way that was less dumb.
        */

        let product_degree = rhs.degree() + self.degree();
        let mut left_matrix_rows: Vec<Vec<Q>> = vec![vec![Q::zero(); self.degree() + 1]; product_degree + 1];

        // Populate left matrix, column by column
        for column_number in 0..(self.degree() + 1) {
            for row_number in 0..(product_degree + 1) {
                if row_number >= column_number && row_number < rhs.degree() + column_number + 1 {
                    left_matrix_rows[row_number][column_number] = rhs.coefficients[row_number - column_number].clone();
                } else {
                    left_matrix_rows[row_number][column_number] = Q::zero();
                }
            }
        }

        // Perform matrix multiplication
        let mut product_coefficients: Vec<Q> = vec![Q::zero(); product_degree + 1];
        for row_number in 0..(product_degree + 1) {
            let mut sum = Q::zero();
            for column_number in 0..(self.degree() + 1) {
                sum += &left_matrix_rows[row_number][column_number] * &self.coefficients[column_number]
            }
            product_coefficients[row_number] = sum;
        }

        Polynomial::truncate_coefficients(&mut product_coefficients);
        self.coefficients = product_coefficients;
    }

    pub fn mul(&self, rhs: &Self) -> Polynomial {
        let mut clone = self.clone();
        clone.mul_mut(rhs);
        clone
    }
}

impl Mul for Polynomial {
    type Output = Polynomial;

    fn mul(self, rhs: Self) -> Self::Output {
        (&self).mul(&rhs)
    }
}

impl Mul<&Self> for Polynomial {
    type Output = Polynomial;

    fn mul(self, rhs: &Self) -> Self::Output {
        (&self).mul(rhs)
    }
}

impl MulAssign for Polynomial {
    fn mul_assign(&mut self, rhs: Self) {
        self.mul_mut(&rhs);
    }
}

impl One for Polynomial {
    fn one() -> Self {
        Polynomial::new(vec![Q::one()])
    }

    fn is_one(&self) -> bool where
        Self: PartialEq, {
        self.eq(&one())
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_multiplication() {
        // t^2 + 2t - 7 * t - 2 == t^3 - 11t + 14
        assert_eq!(Polynomial::from(vec![-2, 1])
                       * Polynomial::from(vec![-7, 2, 1]),
                   Polynomial::from(vec![14, -11, 0, 1]));

        // t^2 - 3t - 10 * t + 2 == t^3 - t^2 - 16t - 20
        assert_eq!(Polynomial::from(vec![-10, -3, 1])
                       * Polynomial::from(vec![2, 1]),
                   Polynomial::from(vec![-20, -16, -1, 1]));

        // 2t^3 - 7t^2 + 4 * t^2 - 1 == 2t^5 - 7t^4 - 2t^3 + 11t^2 - 4
        assert_eq!(Polynomial::from(vec![4, 0, -7, 2])
                       * Polynomial::from(vec![-1, 0, 1]),
                   Polynomial::from(vec![-4, 0, 11, -2, -7, 2]));

        // t^2 + 2t - 7 * t^2 + 2t - 7 == t^4 + 4t^3 - 10t^2 - 28t + 49
        assert_eq!(Polynomial::from(vec![-7, 2, 1])
                       * Polynomial::from(vec![-7, 2, 1]),
                   Polynomial::from(vec![49, -28, -10, 4, 1]));

        // 2t^2 + 2t - 2 * -3t^3 + 2t - 7 == -6t^5 - 6t^4 +10t^3 - 10t^2 - 18t + 14
        assert_eq!(Polynomial::from(vec![-2, 2, 2])
                       * Polynomial::from(vec![-7, 2, 0, -3]),
                   Polynomial::from(vec![14, -18, -10, 10, -6, -6]));
    }

    #[test]
    fn test_multiplication_mut() {
        // t^2 + 2t - 7 * t - 2 == t^3 - 11t + 14
        let mut p1 = Polynomial::from(vec![-2, 1]);
        p1 *= Polynomial::from(vec![-7, 2, 1]);
        assert_eq!(p1, Polynomial::from(vec![14, -11, 0, 1]));

        // t^2 - 3t - 10 * t + 2 == t^3 - t^2 - 16t - 20
        let mut p2 = Polynomial::from(vec![-10, -3, 1]);
        p2 *= Polynomial::from(vec![2, 1]);
        assert_eq!(p2, Polynomial::from(vec![-20, -16, -1, 1]));

        // 2t^3 - 7t^2 + 4 * t^2 - 1 == 2t^5 - 7t^4 - 2t^3 + 11t^2 - 4
        let mut p3 = Polynomial::from(vec![4, 0, -7, 2]);
        p3 *= Polynomial::from(vec![-1, 0, 1]);
        assert_eq!(p3, Polynomial::from(vec![-4, 0, 11, -2, -7, 2]));

        // t^2 + 2t - 7 * t^2 + 2t - 7 == t^4 + 4t^3 - 10t^2 - 28t + 49
        let mut p4 = Polynomial::from(vec![-7, 2, 1]);
        p4 *= Polynomial::from(vec![-7, 2, 1]);
        assert_eq!(p4, Polynomial::from(vec![49, -28, -10, 4, 1]));

        // 2t^2 + 2t - 2 * -3t^3 + 2t - 7 == -6t^5 - 6t^4 +10t^3 - 10t^2 - 18t + 14
        let mut p5 = Polynomial::from(vec![-2, 2, 2]);
        p5 *= Polynomial::from(vec![-7, 2, 0, -3]);
        assert_eq!(p5, Polynomial::from(vec![14, -18, -10, 10, -6, -6]));
    }
}