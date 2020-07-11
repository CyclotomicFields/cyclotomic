use std::ops::{Mul, MulAssign};
use num::{One, one, Zero};

use crate::polynomial::polynomial::{Polynomial, Z, Q};
use std::cmp::max;

impl Polynomial {
    fn coefficient_mul_naive(lhs: &Vec<Q>, rhs: &Vec<Q>) -> Vec<Q> {
        let mut product_coefficients = vec![Q::zero(); lhs.len() + rhs.len() - 1];
        for i in 0..lhs.len() {
            for j in 0..rhs.len() {
                product_coefficients[i + j] += &lhs[i] * &rhs[j];
            }
        }
        product_coefficients
    }

    pub fn mul_mut_naive(&mut self, rhs: &Self) {
        self.coefficients = Polynomial::coefficient_mul_naive(&self.coefficients, &rhs.coefficients);
    }

    fn coefficient_mul_convolutions(lhs: &Vec<Q>, rhs: &Vec<Q>) -> Vec<Q> {
        /*
        Using convolutions

        Algorithm written out clearly here: http://www.cse.ust.hk/~dekai/271/notes/L03/L03.pdf
        */
        let n: usize = max(lhs.len(), rhs.len()) - 1;
        if n < 3 {
            return Polynomial::coefficient_mul_naive(lhs, rhs);
        }
        let m = (n as f64 / 2.0).ceil() as usize;

        #[inline]
        fn split_coefficients(splitting_degree: usize, coefficients: &Vec<Q>) -> (Vec<Q>, Vec<Q>) {
            let (lower, upper) = coefficients.split_at(splitting_degree);
            return (Vec::from(lower), Vec::from(upper));
        }

        let (a0, a1) = split_coefficients(m, &lhs);
        let (b0, b1) = split_coefficients(m, &rhs);

        let mut a0_plus_a1: Vec<Q> = vec![Q::zero(); m + 1];
        for i in 0..m + 1 {
            a0_plus_a1[i] = (a0.get(i).unwrap_or(&Q::zero()) + a1.get(i).unwrap_or(&Q::zero()))
        }
        let mut b0_plus_b1: Vec<Q> = vec![Q::zero(); m + 1];
        for i in 0..m + 1 {
            b0_plus_b1[i] = (b0.get(i).unwrap_or(&Q::zero()) + b1.get(i).unwrap_or(&Q::zero()))
        }
        /*
        Splits up the multiplication over n coefficients into 3 multiplications of n/2 coefficients.
        The naive coefficient-wise multiplication is O(n^2), which is why this is effective.
        */
        let a0b0 = Polynomial::coefficient_mul_convolutions(&a0, &b0);
        let a1b1 = Polynomial::coefficient_mul_convolutions(&a1, &b1);
        let mut a0b1_plus_a1b0 = Polynomial::coefficient_mul_convolutions(&a0_plus_a1, &b0_plus_b1);
        let zero = Q::zero();
        for i in 0..n + 1 {
            a0b1_plus_a1b0[i] -= (a0b0.get(i).unwrap_or(&zero) + a1b1.get(i).unwrap_or(&zero))
        }

        let first_part_coefficients = a0b0;
        let mut mid_part_coefficients = vec![Q::zero(); m];
        mid_part_coefficients.extend(a0b1_plus_a1b0);
        let mut end_part_coefficients = vec![Q::zero(); 2 * m];
        end_part_coefficients.extend(a1b1);
        for i in 0..end_part_coefficients.len() {
            end_part_coefficients[i] += (first_part_coefficients.get(i).unwrap_or(&Q::zero()) + mid_part_coefficients.get(i).unwrap_or(&Q::zero()));
        }
        end_part_coefficients
    }

    pub fn mul_mut_convolutions(&mut self, rhs: &Self) {
        self.coefficients = Polynomial::coefficient_mul_convolutions(&self.coefficients, &rhs.coefficients);
    }

    pub fn mul_mut_fft(&mut self, _rhs: &Self) {
        /*
        Using the Fast Fourier Transform

        Algorithm written out clearly here: https://agarri.ga/post/multiply-polynomials-n-log-n-time/
        */
    }

    pub fn mul(&self, rhs: &Self) -> Polynomial {
        let mut clone = self.clone();
        clone.mul_mut_convolutions(rhs);
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
        self.mul_mut_convolutions(&rhs);
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

        // 2 + 5t + 3t^2 + t^3 - t^4 * 1 + 2t + 2t^2 + 3t^3 + 6t^4 == 2 + 9t + 17t^2 + 23t^3 + 34t^4 + 39t^5 + 19t^6 + 3t^7 - 6t^8
        assert_eq!(Polynomial::from(vec![2, 5, 3, 1, -1])
                       * Polynomial::from(vec![1, 2, 2, 3, 6]),
                   Polynomial::from(vec![2, 9, 17, 23, 34, 39, 19, 3, -6]));
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

        // 2 + 5t + 3t^2 + t^3 - t^4 * 1 + 2t + 2t^2 + 3t^3 + 6t^4 == 2 + 9t + 17t^2 + 23t^3 + 34t^4 + 39t^5 + 19t^6 + 3t^7 - 6t^8
        let mut p6 = Polynomial::from(vec![2, 5, 3, 1, -1]);
        p6 *= Polynomial::from(vec![1, 2, 2, 3, 6]);
        assert_eq!(p6, Polynomial::from(vec![2, 9, 17, 23, 34, 39, 19, 3, -6]));
    }
}