use std::ops::Mul;

use num::{One, one, Zero};

use crate::polynomial::polynomial::{Polynomial, Z};

impl Polynomial {
    // TODO: Figure out how to use a fast library for this.
    pub fn mul(&self, rhs: &Self) -> Polynomial {
        /*
        Matrix product method

        Represent the coefficients of the lhs polynomial as diagonal entries in
        a not-necessarily-square matrix, and perform a matrix multiplication
        with a column vector containing the coefficients of the rhs polynomial.
        */

        let lhs = self;
        let product_degree = lhs.degree() + rhs.degree();
        let mut left_matrix_rows: Vec<Vec<Z>> = vec![vec![Z::zero(); rhs.degree() + 1]; product_degree + 1];

        // Populate left matrix, column by column
        for column_number in 0..(rhs.degree() + 1) {
            for row_number in 0..(product_degree + 1) {
                if row_number >= column_number && row_number < lhs.degree() + column_number + 1 {
                    left_matrix_rows[row_number][column_number] = lhs.coefficients[row_number - column_number].clone();
                } else {
                    left_matrix_rows[row_number][column_number] = Z::zero();
                }
            }
        }

        // Perform matrix multiplication
        let mut product_coefficients: Vec<Z> = vec![Z::zero(); product_degree + 1];
        for row_number in 0..(product_degree + 1) {
            let mut sum = Z::zero();
            for column_number in 0..(rhs.degree() + 1) {
                sum += &left_matrix_rows[row_number][column_number] * &rhs.coefficients[column_number]
            }
            product_coefficients[row_number] = sum;
        }

        Polynomial::truncate_coefficients(&mut product_coefficients);
        Polynomial::new(product_coefficients)
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

impl One for Polynomial {
    fn one() -> Self {
        Polynomial::new(vec![Z::one()])
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
}