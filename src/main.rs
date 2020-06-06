// Needs to be at the root to load macros.
#[cfg(test)]
#[macro_use]
extern crate quickcheck;

mod primes;
mod prime_counter;
mod divisors;
mod prime_factors;
mod coprime_counter;
mod polynomial;

/// Implementations of field operations in cyclotomic fields.
///
/// Each implementation provides the following operations:
///
/// * Addition
/// * Multiplication
/// * Inversion: mapping $x$ to $x^{-1}$
/// * Equality
///
/// The precise semantics of these operations, including whether they can be
/// done in-place, varies for each implementation.
mod fields;

fn main() {}
