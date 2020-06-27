use std::ops::{Add, Neg, Sub};

use crate::polynomial::polynomial::Polynomial;

impl Polynomial {
    // TODO: Stop griefing
    pub fn sub(&self, rhs: &Self) -> Polynomial {
        self.clone().add(rhs.neg())
    }
}

impl Sub for Polynomial {
    type Output = Polynomial;

    fn sub(self, rhs: Self) -> Self::Output {
        (&self).sub(&rhs)
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_subtraction() {
        // t^2 - 3t - 10 - (t + 2) == t^2 - 4t - 12
        assert_eq!(Polynomial::from(vec![-10, -3, 1])
                       - Polynomial::from(vec![2, 1]),
                   Polynomial::from(vec![-12, -4, 1]));

        // t^2 + 2t - 7 - (t - 2) == t^2 + t - 5
        assert_eq!(Polynomial::from(vec![-7, 2, 1])
                       - Polynomial::from(vec![-2, 1]),
                   Polynomial::from(vec![-5, 1, 1]));

        // t^2 + 2t - 7 - (t^2 + 2t - 7) == 0
        assert_eq!(Polynomial::from(vec![-7, 2, 1])
                       - Polynomial::from(vec![-7, 2, 1]),
                   Polynomial::from(vec![0]));

        // -3t^3 + 2t - 7 - (2t^2 + 2t - 2) == -3t^3 - 2t^2 - 5
        assert_eq!(Polynomial::from(vec![-7, 2, 0, -3])
                       - Polynomial::from(vec![-2, 2, 2]),
                   Polynomial::from(vec![-5, 0, -2, -3]));
    }
}