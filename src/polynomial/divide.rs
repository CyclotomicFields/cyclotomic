use std::ops::{Div, MulAssign};

use num::{Integer, One, Zero};

use crate::polynomial::polynomial::{Polynomial, Q, Z};

impl Polynomial {
    pub fn div(&self, divisor: &Self) -> (Polynomial, Polynomial) {
        /*
        Polynomial Long Division

        We have a dividend and a divisor, which are both polynomials, where the
        divisor has degree less than or equal to the dividend. These are both
        non-constant polynomials, that is, they have a degree of at least 1.

        We are seeking a quotient and a remainder for the division of these two
        polynomials. These are both polynomials, where the quotient has degree
        strictly less than the degree of the dividend, and the remainder has
        degree strictly less than the divisor.

        The algorithm involves iterating on the dividend. We shall call the
        polynomial that develops through these iterations as the current
        polynomial. In each iteration, we add terms to the polynomial that will
        end up by the quotient. We shall refer to this in-progress polynomial
        as the quotient accumulator.

        For each iteration, perform the below steps.

        First, divide the leading term of the current polynomial by the
        leading term of the divisor. Initially, as mentioned earlier, the
        current polynomial is the dividend. This result we'll' call the
        current term. Append the current term to the quotient accumulator.

        We now perform steps to iterate on the current polynomial. Multiply the
        current term by the divisor, and the subtract that from the current
        polynomial to get the next iteration. Below is some pseudocode to
        hopefully make this clearer.

        current polynomial -= current term * divisor

        We continue iterating until the the degree of the remainder is strictly
        less than the degree of the divisor.
        */
        assert!(self.degree() >= divisor.degree());
        assert!(!divisor.is_zero());
        if divisor.is_one() {
            return (self.clone(), Polynomial::zero());
        } else if divisor.neg().is_one() {
            return (self.neg(), Polynomial::zero());
        } else if self.degree() == 0 && divisor.degree() == 0 {
            let rational_number = self.leading_term_coefficient().div(divisor.leading_term_coefficient());
            return (Polynomial::new(vec![rational_number]), Polynomial::zero());
        } else if divisor.degree() == 0 {
            let divided_coefficients = self.coefficients.iter().map(|c| c.div(&divisor.coefficients[0])).collect();
            return (Polynomial::new(divided_coefficients), Polynomial::zero());
        }
        assert!(self.degree() > 0);
        let mut quotient_accumulator = vec![Q::zero(); self.degree()];
        let mut current_polynomial = self.clone();

        while current_polynomial.degree() >= divisor.degree() {
            let current_term_degree = current_polynomial.degree() - divisor.degree();
            let current_term_coefficient = current_polynomial.leading_term_coefficient().div(divisor.leading_term_coefficient());
            quotient_accumulator[current_term_degree] = current_term_coefficient.clone();
            let mut subtraction_coefficients = vec![Q::zero(); current_term_degree];
            subtraction_coefficients.extend(divisor.coefficients.iter().map(|c| Q::from(c.clone())));
            subtraction_coefficients.iter_mut().for_each(|c| c.mul_assign(current_term_coefficient.clone()));
            current_polynomial = current_polynomial - Polynomial::new(subtraction_coefficients);
        }

        Polynomial::truncate_coefficients(&mut quotient_accumulator);
        (Polynomial::new(quotient_accumulator), current_polynomial)
    }
}

impl Div<&Self> for Polynomial {
    type Output = (Polynomial, Polynomial);

    fn div(self, divisor: &Self) -> Self::Output {
        (&self).div(divisor)
    }
}

#[cfg(test)]
mod polynomial_tests {
    use num::{One, Zero};
    use super::*;

    #[inline]
    fn vec_q(numerators: Vec<i64>, denominators: Vec<i64>) -> Vec<Q> {
        let mut qs = vec![];
        for i in 0..numerators.len() {
            qs.push(Q::new(Z::from(numerators[i]), Z::from(denominators[i])));
        }
        qs
    }

    #[test]
    fn test_long_division() {
        // t^2 + 6 / -6t - 1 == -(1/6)t + (1/36) remainder (217/36)
        assert_eq!(Polynomial::from(vec![6, 0, 1])
                       .div(&Polynomial::from(vec![-1, -6])),
                   (Polynomial::new(vec_q(vec![1, -1], vec![36, 6])),
                    Polynomial::new(vec_q(vec![217], vec![36]))));

        // 14 / 6 == 7 / 3 remainder 0
        assert_eq!(Polynomial::from(vec![14])
                       .div(&Polynomial::from(vec![6])),
                   (Polynomial::new(vec_q(vec![7], vec![3])), Polynomial::zero()));

        // t^2 - 1 / 1 == t^2 - 1 remainder 0
        assert_eq!(Polynomial::from(vec![-1, 0, 1])
                       .div(&Polynomial::from(vec![1])),
                   (Polynomial::from(vec![-1, 0, 1]), Polynomial::zero()));

        // t + 2 / -1 == -t - 2 remainder 0
        assert_eq!(Polynomial::from(vec![2, 1])
                       .div(&Polynomial::from(vec![-1])),
                   (Polynomial::from(vec![-2, -1]), Polynomial::zero()));

        // t^2 - 3t - 10 / t + 2 == t - 5 remainder 0
        assert_eq!(Polynomial::from(vec![-10, -3, 1])
                       .div(&Polynomial::from(vec![2, 1])),
                   (Polynomial::from(vec![-5, 1]), Polynomial::zero()));

        // t^2 + 2t - 7 / t - 2 == t + 4 remainder 1
        assert_eq!(Polynomial::from(vec![-7, 2, 1])
                       .div(&Polynomial::from(vec![-2, 1])),
                   (Polynomial::from(vec![4, 1]), Polynomial::one()));

        // 2t^3 - 7t^2 + 4 / t^2 - 1 == 2t - 7 remainder 2t - 3
        assert_eq!(Polynomial::from(vec![4, 0, -7, 2])
                       .div(&Polynomial::from(vec![-1, 0, 1])),
                   (Polynomial::from(vec![-7, 2]), Polynomial::from(vec![-3, 2])));

        // t + 2 / 1 == t + 2 remainder 0
        assert_eq!(Polynomial::from(vec![2, 1])
                       .div(&Polynomial::one()),
                   (Polynomial::from(vec![2, 1]), Polynomial::zero()));
    }
}