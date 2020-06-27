use crate::polynomial::polynomial::{Polynomial, Z, Q};
use std::fmt::{Display, Formatter};
use std::fmt;
use num::{Zero, One, ToPrimitive};
use std::ops::Neg;

impl Polynomial {
    pub fn to_string(&self) -> String {
        if self.is_zero() {
            return format!("{}", 0);
        } else if self.is_one() {
            return format!("{}", 1);
        }

        #[inline]
        fn format_coefficient(coefficient: Q, degree: usize, i: usize) -> String {
            if i == 0 {
                return if coefficient.denom().is_one() {
                    if let Some(c) = coefficient.to_integer().to_i64() {
                        if coefficient.is_one() && degree > 0 {
                            "".to_string()
                        } else if coefficient.eq(&Q::one().neg()) && degree > 0 {
                            "-".to_string()
                        } else {
                            format!("{}", c)
                        }
                    } else {
                        format!("({})", coefficient.to_string())
                    }
                } else {
                    format!("({})", coefficient.to_string())
                };
            }
            if coefficient.is_one() && degree > 0 {
                "+ ".to_string()
            } else if coefficient.eq(&Q::one().neg()) && degree > 0 {
                "- ".to_string()
            } else {
                return if coefficient.denom().is_one() {
                    if let Some(c) = coefficient.to_integer().to_i64() {
                        if c < 0 {
                            format!("- {}", c * -1)
                        } else {
                            format!("+ {}", c)
                        }
                    } else {
                        format!("({})", coefficient.to_string())
                    }
                } else {
                    format!("({})", coefficient.to_string())
                }
            }
        }

        #[inline]
        fn format_variable(degree: usize) -> String {
            if degree == 0 {
                "".to_string()
            } else if degree == 1 {
                "t".to_string()
            } else {
                format!("t^{}", degree)
            }
        }

        let mut string_builder = "".to_string();
        for i in 0..self.coefficients.len() {
            let coefficient = self.coefficients[self.coefficients.len() - 1 - i].clone();
            if !coefficient.is_zero() {
                let degree = self.degree() - i;
                string_builder.push_str(
                    format!("{}{} ",
                            format_coefficient(coefficient, degree, i),
                            format_variable(degree))
                        .as_str())
            }
        }
        string_builder.truncate(string_builder.len() - 1);
        string_builder
    }
}

impl Display for Polynomial {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        return write!(f, "{}", self.to_string());
    }
}

#[cfg(test)]
mod polynomial_tests {
    use super::*;

    #[test]
    fn test_to_string() {
        assert_eq!(Polynomial::from(vec![-1]).to_string(), "-1".to_string());
        assert_eq!(Polynomial::one().to_string(), "1".to_string());
        assert_eq!(Polynomial::zero().to_string(), "0".to_string());
        assert_eq!(Polynomial::from(vec![-2, 1]).to_string(), "t - 2".to_string());
        assert_eq!(Polynomial::from(vec![1, 0, -7, 2]).to_string(), "2t^3 - 7t^2 + 1".to_string());
        assert_eq!(Polynomial::from(vec![-1, 2]).to_string(), "2t - 1".to_string());
        assert_eq!(Polynomial::from(vec![0, -1]).to_string(), "-t".to_string());
        assert_eq!(Polynomial::from(vec![2, 5, 4, 6]).to_string(), "6t^3 + 4t^2 + 5t + 2".to_string());
        assert_eq!(Polynomial::from(vec![2, 5, 0, 0, -9]).to_string(), "-9t^4 + 5t + 2".to_string());
        assert_eq!(Polynomial::from(vec![1, 1, 1]).to_string(), "t^2 + t + 1".to_string());
        assert_eq!(Polynomial::from(vec![1, -1, 1]).to_string(), "t^2 - t + 1".to_string());
    }
}