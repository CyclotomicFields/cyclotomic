type R = f64;
type ZPlus = u64;

pub trait PrimeCounter {
    fn pi(&self, x: R) -> ZPlus;
}