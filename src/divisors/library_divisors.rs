extern crate num;
extern crate divisors;

use self::divisors::get_divisors;
use self::num::{BigInt, BigRational, Zero, ToPrimitive};
use crate::divisors::divisors::Divisors;

type Z = num::bigint::BigInt;

pub struct LibraryDivisors {}

impl LibraryDivisors {
    pub fn new() -> LibraryDivisors {
        LibraryDivisors {}
    }
}

impl Divisors for LibraryDivisors {
    fn divisors(&self, z: Z) -> Vec<Z> {
        if z.is_zero() {
            panic!("Can't get the divisors of 0")
        }
        let mut z_i128 = z.to_i128().unwrap();
        if z_i128 < 1 {
            z_i128 = -1 * z_i128
        }
        let z_u128 = z_i128 as u128;
        let mut divisors = vec![1];
        divisors.extend(get_divisors(z_u128));
        if (z_u128 != 1) && (z_u128 != 2) {
            divisors.push(z_u128);
        }
        return divisors.iter().map(|&u| Z::from(u)).collect();
    }

    fn divisors_without_one(&self, z: Z) -> Vec<Z> {
        if z.is_zero() {
            panic!("Can't get the divisors of 0")
        }
        let mut z_i128 = z.to_i128().unwrap();
        if z_i128 < 1 {
            z_i128 = -1 * z_i128
        }
        let z_u128 = z_i128 as u128;
        let mut divisors = vec![];
        divisors.extend(get_divisors(z_u128));
        if (z_u128 != 1) && (z_u128 != 2) {
            divisors.push(z_u128);
        }
        return divisors.iter().map(|&u| Z::from(u)).collect();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_divisors() {
        let strategy = LibraryDivisors::new();
        assert_eq!(strategy.divisors(Z::from(2)), vec_z(vec![1, 2]));
        assert_eq!(strategy.divisors(Z::from(6)), vec_z(vec![1, 2, 3, 6]));
        assert_eq!(strategy.divisors(Z::from(1)), vec_z(vec![1]));
        assert_eq!(strategy.divisors(Z::from(12)), vec_z(vec![1, 2, 3, 4, 6, 12]));
        assert_eq!(strategy.divisors(Z::from(-10)), vec_z(vec![1, 2, 5, 10]));
    }

    #[test]
    fn test_divisors_without_one() {
        let strategy = LibraryDivisors::new();
        assert_eq!(strategy.divisors_without_one(Z::from(2)), vec_z(vec![2]));
        assert_eq!(strategy.divisors_without_one(Z::from(4)), vec_z(vec![2, 4]));
        assert_eq!(strategy.divisors_without_one(Z::from(6)), vec_z(vec![2, 3, 6]));
        assert_eq!(strategy.divisors_without_one(Z::from(1)), vec_z(vec![]));
        assert_eq!(strategy.divisors_without_one(Z::from(12)), vec_z(vec![2, 3, 4, 6, 12]));
        assert_eq!(strategy.divisors_without_one(Z::from(-10)), vec_z(vec![2, 5, 10]));
    }

    fn vec_z(vec: Vec<i64>) -> Vec<Z> {
        vec.iter().map(|&i| BigInt::from(i)).collect()
    }
}