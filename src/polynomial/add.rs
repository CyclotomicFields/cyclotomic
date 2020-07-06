use std::ops::{Add, AddAssign};

use num::{Zero, zero};

use crate::polynomial::polynomial::{Polynomial, Z, Q};

impl Add for Polynomial {
    type Output = Polynomial;

    fn add(self, rhs: Self) -> Self::Output {
        // Ensure that self's degree is greater than or equal to rhs's degree
        if rhs.degree() > self.degree() {
            return rhs.add(self);
        }

        let mut new_coefficients = self.coefficients.clone();
        let mut i = 0;
        while i < rhs.coefficients.len() {
            new_coefficients[i] = &self.coefficients[i] + &rhs.coefficients[i];
            i += 1;
        }
        while i < self.coefficients.len() {
            new_coefficients[i] = self.coefficients[i].clone();
            i += 1;
        }
        Polynomial::truncate_coefficients(&mut new_coefficients);
        Polynomial::new(new_coefficients)
    }
}

impl Zero for Polynomial {
    fn zero() -> Self {
        Polynomial::new(vec![Q::zero()])
    }

    fn is_zero(&self) -> bool {
        self.eq(&zero())
    }
}

impl AddAssign for Polynomial {
    fn add_assign(&mut self, rhs: Self) {
        unimplemented!()
    }
}

// TODO: Use quickcheck instead
#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_addition() {
        // t^2 + 2t - 7 (+) t - 2
        assert_eq!(Polynomial::from(vec![-2, 1])
                       .add(Polynomial::from(vec![-7, 2, 1])),
                   Polynomial::from(vec![-9, 3, 1]));

        // t^2 - t - 10 (+) t + 2
        assert_eq!(Polynomial::from(vec![-10, -1, 1])
                       .add(Polynomial::from(vec![2, 1])),
                   Polynomial::from(vec![-8, 0, 1]));

        // 2t^3 - 7t^2 + 1 (+) t^2 - 1
        assert_eq!(Polynomial::from(vec![1, 0, -7, 2])
                       .add(Polynomial::from(vec![-1, 0, 1])),
                   Polynomial::from(vec![0, 0, -6, 2]));

        // t^2 + 2t - 7 (+) -t^2 + 2t - 7
        assert_eq!(Polynomial::from(vec![-7, 2, 1])
                       .add(Polynomial::from(vec![-7, 2, -1])),
                   Polynomial::from(vec![-14, 4]));
    }

    #[test]
    #[ignore]
    fn test_mut_addition() {
        // t^2 + 2t - 7 (+) t - 2
        let mut p1 = Polynomial::from(vec![-2, 1]);
        p1.add_assign(Polynomial::from(vec![-7, 2, 1]));
        assert_eq!(p1, Polynomial::from(vec![-9, 3, 1]));

        // t^2 - t - 10 (+) t + 2
        let mut p2 = Polynomial::from(vec![-10, -1, 1]);
        p2.add_assign(Polynomial::from(vec![2, 1]));
        assert_eq!(p2, Polynomial::from(vec![-8, 0, 1]));

        // 2t^3 - 7t^2 + 1 (+) t^2 - 1
        let mut p3 = Polynomial::from(vec![1, 0, -7, 2]);
        p3.add_assign(Polynomial::from(vec![-1, 0, 1]));
        assert_eq!(p3, Polynomial::from(vec![0, 0, -6, 2]));

        // t^2 + 2t - 7 (+) -t^2 + 2t - 7
        let mut p4 = Polynomial::from(vec![-7, 2, 1]);
        p4.add_assign(Polynomial::from(vec![-7, 2, -1]));
        assert_eq!(p4, Polynomial::from(vec![-14, 4]));
    }
}