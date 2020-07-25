type R = f64;
type ZPlus = u64;

/*
A trait for objects which can perform pi(x), which counts the number of primes
not larger than x. This function is needed for fast algorithms for factorizing
numbers.
*/
pub trait PrimeCounter {
    fn pi(&self, x: R) -> ZPlus;
}
