use crate::fields::Z;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};
use rug::ops::Pow;

pub trait Power {
    fn power(&self, i: u32) -> Self;
}

/// Something suitable to be the i in e^i. Big or small signed ints, most likely,
/// and some other useful numeric functions.
pub trait Exponent:
    Clone
    + Add<Output = Self>
    + Mul<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + Neg<Output = Self>
    + From<i64>
    + Eq
    + Display
    + Send
    + Hash
    + Debug
    + Power
    + Sync
{
    fn gcd(x: &Self, y: &Self) -> Self;
    fn lcm(x: &Self, y: &Self) -> Self;

    /// This doesn't support some numbers bigger than 2^(2^32).
    fn factorise(n: &Self) -> Vec<(Self, u32)> {
        let mut result = vec![];
        let mut n_factored = n.clone();
        let mut current_divisor = Self::from(2);

        while current_divisor != *n {
            let mut power: u32 = 0;

            while (n_factored.clone() % current_divisor.clone()) == Self::from(0) {
                power += 1;
                n_factored = n_factored / current_divisor.clone();
            }

            if power != 0 {
                result.push((current_divisor.clone(), power));
            }

            // TODO: this is really bad, improve it later
            current_divisor = current_divisor + Self::from(1);
        }

        if result.is_empty() {
            // This means n has no divisors smaller than itself other than 1,
            // so is prime.
            result.push((n.clone(), 1))
        }

        result
    }

    fn phi(n: &Self) -> Self {
        let mut count = Self::from(0);
        let mut k = Self::from(1);
        while &k != n {
            if Exponent::gcd(n, &k) == Self::from(1) {
                count = count + Self::from(1);
            }
            k = k + Self::from(1);
        }
        count
    }

    fn math_mod(x: &Self, n: &Self) -> Self {
        (x.clone() % n.clone() + n.clone()) % n.clone()
    }
}

impl Exponent for Z {
    fn gcd(x: &Z, y: &Z) -> Z {
        x.clone().gcd(y)
    }

    fn lcm(x: &Z, y: &Z) -> Z {
        x.clone().lcm(y)
    }
}

impl Power for Z {
    fn power(&self, i: u32) -> Self {
        self.pow(i).into()
    }
}

impl Exponent for i64 {
    fn gcd(x: &i64, y: &i64) -> i64 {
        num::integer::gcd(*x as u64, *y as u64) as i64
    }

    // There is a real danger that the lcm doesn't fit into
    // a fixed-size integer, but we ignore that.
    fn lcm(x: &i64, y: &i64) -> i64 {
        num::integer::lcm(*x as u64, *y as u64) as i64
    }
}

impl Power for i64 {
    fn power(&self, i: u32) -> Self {
        self.pow(i)
    }
}