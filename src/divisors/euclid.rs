extern crate num;

use self::num::{BigInt, BigRational, Zero, ToPrimitive, One, Integer};
use std::ops::{Div, Rem};

type Z = num::bigint::BigInt;

pub struct Euclid {}

impl Euclid {
    pub fn new() -> Euclid {
        Euclid {}
    }

    pub fn gcd(&self, m: Z, n: Z) -> Z {
        if m == n {
            return m;
        } else if n > m {
            return self.gcd(n, m);
        }
        let mut trailing_row = vec![m, Z::one(), Z::zero()];
        let mut leading_row = vec![n, Z::zero(), Z::one()];
        while !&leading_row[0].is_zero() {
            let (quotient, remainder) = (&trailing_row[0]).div_rem(&leading_row[0]);
            let new_leading_row = vec![remainder,
                                       &trailing_row[1] - &quotient * &leading_row[1],
                                       &trailing_row[2] - &quotient * &leading_row[2]];
            trailing_row = leading_row;
            leading_row = new_leading_row;
        }
        return trailing_row[0].clone();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extended_euclid_gcd() {
        let euclid = Euclid::new();
        assert_eq!(euclid.gcd(Z::from(80), Z::from(62)), Z::from(2));
        assert_eq!(euclid.gcd(Z::from(5), Z::from(5)), Z::from(5));
        assert_eq!(euclid.gcd(Z::from(420), Z::from(1782)), Z::from(6));
        assert_eq!(euclid.gcd(Z::from(371250), Z::from(873630)), Z::from(90));
        assert_eq!(euclid.gcd(Z::from(3283741336253_i64), Z::from(0234551825263_i64)), Z::from(7));
    }
}