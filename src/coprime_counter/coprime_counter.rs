type ZPlus = u64;

/*
A trait for structs which can compute phi(n), which counts the integers that
are coprime with n. This is known as Euler's totient function, denoted with a
phi.
*/
pub trait CoprimeCounter {
    fn phi(&self, n: ZPlus) -> ZPlus;
}