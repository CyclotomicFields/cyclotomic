extern crate num;

use num::pow::pow;

type Z = i128;
type ZPlus = usize;

pub struct Polynomial {
    coefficients: Vec<Z>,
    degrees: Vec<ZPlus>,
}

impl Polynomial {
    pub fn new(coefficients: Vec<Z>, degrees: Vec<ZPlus>) -> Polynomial {
        assert_eq!(coefficients.len(), degrees.len());
        Polynomial { coefficients, degrees }
    }

    pub fn substitute(&self, t: Z) -> Z {
        let mut sum: Z = 0;
        for j in 0..self.coefficients.len() {
            sum += self.coefficients[j] * pow(t, self.degrees[j]);
        }
        return sum;
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_polynomial_substitution() {
        // p = 1
        let mut p: Polynomial = Polynomial::new(vec![1], vec![0]);

        assert_eq!(p.substitute(5), 1);
        assert_eq!(p.substitute(0), 1);
        assert_eq!(p.substitute(-1), 1);
        assert_eq!(p.substitute(-14), 1);
        assert_eq!(p.substitute(1), 1);
        assert_eq!(p.substitute(2), 1);
        assert_eq!(p.substitute(125716), 1);

        // p = 2t + 1
        p = Polynomial::new(vec![1, 2], vec![0, 1]);

        assert_eq!(p.substitute(5), 11);
        assert_eq!(p.substitute(0), 1);
        assert_eq!(p.substitute(-1), -1);
        assert_eq!(p.substitute(-14), -27);
        assert_eq!(p.substitute(1), 3);
        assert_eq!(p.substitute(2), 5);
        assert_eq!(p.substitute(125716), 251433);

        // p = t^5 - t^3 - 2t^2 - 1
        p = Polynomial::new(vec![1, 0, -1, -2, 0, -1], vec![5, 4, 3, 2, 1, 0]);

        assert_eq!(p.substitute(5), 2949);
        assert_eq!(p.substitute(0), -1);
        assert_eq!(p.substitute(-1), -3);
        assert_eq!(p.substitute(-14), -535473);
        assert_eq!(p.substitute(1), -3);
        assert_eq!(p.substitute(2), 15);
        assert_eq!(p.substitute(125716), 31401671890851373618737567);
    }
}