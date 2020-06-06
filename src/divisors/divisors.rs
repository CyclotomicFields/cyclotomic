extern crate num;

use self::num::{BigInt};

type Z = num::bigint::BigInt;

pub trait Divisors {
    /*
    Returns a vector containing all the divisors of z, including 1 and itself.
    */
    fn divisors(&self, z: &Z) -> Vec<Z>;

    /*
    Returns a vector containing all the divisors of z, including itself, but not 1.
    */
    fn divisors_without_one(&self, z: &Z) -> Vec<Z>;
}