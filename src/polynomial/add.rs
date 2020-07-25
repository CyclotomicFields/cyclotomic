use std::ops::{Add, AddAssign};

use num::{zero, Zero};

use crate::polynomial::polynomial::{Polynomial, Q};

impl Polynomial {
    fn add_mut(&mut self, rhs: Self) {
        if self.degree() >= rhs.degree() {
            let mut i = 0;
            while i < rhs.coefficients.len() {
                self.coefficients[i] += &rhs.coefficients[i];
                i += 1;
            }
        } else {
            let mut i = 0;
            while i < self.coefficients.len() {
                self.coefficients[i] += &rhs.coefficients[i];
                i += 1;
            }
            // This is suboptimal. It could be implemented faster by doing a
            // bulk copy of the slice of data corresponding to the additional
            // coefficients present in the RHS polynomial.
            self.coefficients
                .extend(vec![Q::zero(); rhs.degree() - self.degree()]);
            while i < rhs.coefficients.len() {
                self.coefficients[i] = rhs.coefficients[i].clone();
                i += 1;
            }
        }
        Polynomial::truncate_coefficients(&mut self.coefficients)
    }
}

impl Add for Polynomial {
    type Output = Polynomial;

    fn add(self, rhs: Self) -> Self::Output {
        // Possibly suboptimal
        let mut clone = self.clone();
        clone.add_mut(rhs);
        clone
    }
}

impl AddAssign for Polynomial {
    fn add_assign(&mut self, rhs: Self) {
        self.add_mut(rhs);
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

// TODO: Use quickcheck instead
#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_addition() {
        // t^2 + 2t - 7 (+) t - 2
        assert_eq!(
            Polynomial::from(vec![-2, 1]).add(Polynomial::from(vec![-7, 2, 1])),
            Polynomial::from(vec![-9, 3, 1])
        );

        // t^2 - t - 10 (+) t + 2
        assert_eq!(
            Polynomial::from(vec![-10, -1, 1]).add(Polynomial::from(vec![2, 1])),
            Polynomial::from(vec![-8, 0, 1])
        );

        // 2t^3 - 7t^2 + 1 (+) t^2 - 1
        assert_eq!(
            Polynomial::from(vec![1, 0, -7, 2]).add(Polynomial::from(vec![-1, 0, 1])),
            Polynomial::from(vec![0, 0, -6, 2])
        );

        // t^2 + 2t - 7 (+) -t^2 + 2t - 7
        assert_eq!(
            Polynomial::from(vec![-7, 2, 1]).add(Polynomial::from(vec![-7, 2, -1])),
            Polynomial::from(vec![-14, 4])
        );
    }

    #[test]
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

        // t^2 - 1 (+) 2t^3 - 7t^2 + 1
        let mut p5 = Polynomial::from(vec![-1, 0, 1]);
        p5.add_assign(Polynomial::from(vec![1, 0, -7, 2]));
        assert_eq!(p5, Polynomial::from(vec![0, 0, -6, 2]));
    }
}
