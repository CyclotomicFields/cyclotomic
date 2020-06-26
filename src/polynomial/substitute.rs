use std::ops::Mul;

use num::Zero;

use crate::polynomial::polynomial::{Polynomial, Q, Z};

impl Polynomial {
    pub fn substitute(&self, t: Q) -> Q {
        let mut sum: Q = Q::zero();
        for j in 0..self.coefficients.len() {
            sum += num::pow(t.clone(), j).mul(&self.coefficients[j]);
        }
        return sum;
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_polynomial_substitution() {
        #[inline]
        fn q_from_i64(n: i64) -> Q {
            return Q::from(Z::from(n));
        }

        // p = 1
        let mut p: Polynomial = Polynomial::from(vec![1]);

        assert_eq!(p.substitute(q_from_i64(5)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(0)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(-1)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(-14)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(1)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(2)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(125716)), q_from_i64(1));

        // p = 2t + 1
        p = Polynomial::from(vec![1, 2]);

        assert_eq!(p.substitute(q_from_i64(5)), q_from_i64(11));
        assert_eq!(p.substitute(q_from_i64(0)), q_from_i64(1));
        assert_eq!(p.substitute(q_from_i64(-1)), q_from_i64(-1));
        assert_eq!(p.substitute(q_from_i64(-14)), q_from_i64(-27));
        assert_eq!(p.substitute(q_from_i64(1)), q_from_i64(3));
        assert_eq!(p.substitute(q_from_i64(2)), q_from_i64(5));
        assert_eq!(p.substitute(q_from_i64(125716)), q_from_i64(251433));

        // p = t^5 - t^3 - 2t^2 - 1
        p = Polynomial::from(vec![-1, 0, -2, -1, 0, 1]);

        assert_eq!(p.substitute(q_from_i64(5)), q_from_i64(2949));
        assert_eq!(p.substitute(q_from_i64(0)), q_from_i64(-1));
        assert_eq!(p.substitute(q_from_i64(-1)), q_from_i64(-3));
        assert_eq!(p.substitute(q_from_i64(-14)), q_from_i64(-535473));
        assert_eq!(p.substitute(q_from_i64(1)), q_from_i64(-3));
        assert_eq!(p.substitute(q_from_i64(2)), q_from_i64(15));
        assert_eq!(p.substitute(q_from_i64(123)), q_from_i64(28151165717_i64));
    }
}