type ZPlus = u64;

/*
A trait for structs which can compute the number of prime factors of an integer
n. Returns them as a vector which may have duplicated values.
*/
pub trait PrimeFactorize {
    fn prime_factors(&self, n: ZPlus) -> Vec<ZPlus>;
}