#![feature(test)]

extern crate rand;
extern crate test;

use cyclotomic::fields::sparse::{random_cyclotomic, Number};
use cyclotomic::fields::AdditiveGroup;
use cyclotomic::fields::FieldElement;
use cyclotomic::fields::MultiplicativeGroup;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use test::{black_box, Bencher};
use std::time::Instant;

fn main() {
    let gen = &mut ChaCha20Rng::seed_from_u64(12345);
    let num_bench = 10000;

    let mut nums: Vec<Number> = (1..num_bench * 6 + 1)
        .into_iter()
        .map(|_| random_cyclotomic(gen))
        .collect();
    let mut i = 0;

    let start = Instant::now();
    for _ in 1..=num_bench {
        let x = &mut nums[i].clone();
        let y = &mut nums[i + 1].clone();
        let z = &mut nums[i + 2].clone();
        let a = &mut nums[i + 3].clone();
        let b = &mut nums[i + 4].clone();
        let c = &mut nums[i + 5].clone();

        // black_box means don't optimise this away into a no-op
        black_box(x.mul(y.add(z)).eq(a.mul(b.add(c))));

        i = i + 6;
    }

    println!("time elapsed: {} ms", start.elapsed().as_millis());
}
