// Based on the Quadratic Sieve algorithm, which is apparently fast: https://en.wikipedia.org/wiki/Quadratic_sieve

use std::ops::Div;

type ZPlus = u64;

// The first step of the algorithm is to choose a smoothness bound. In other words, we are looking
// to decide on a number B, where we believe that no prime factors of n are greater than B.
fn choose_smoothness_bound(n: ZPlus) -> ZPlus {
    return (n as f64).sqrt().ceil() as ZPlus;
}

#[cfg(test)]
mod quadratic_sieve_tests {
    use super::*;

    #[test]
    fn test_choose_smoothness_bound() {
        assert_gte(choose_smoothness_bound(1620), 5, 41);
        assert_gte(choose_smoothness_bound(49), 7, 7);
        assert_gte(choose_smoothness_bound(15750), 7, 126);
        assert_gte(choose_smoothness_bound(702), 13, 27);
        assert_gte(choose_smoothness_bound(121), 11, 11);
    }

    fn assert_gte(a: u64, lower_bound: u64, upper_bound: u64) {
        assert!(a >= lower_bound, "{} was unexpectedly not greater than or equal to {}", a, lower_bound);
        assert!(a <= upper_bound, "{} was unexpectedly not less than or equal to {}", a, upper_bound);
    }
}