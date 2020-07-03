#![feature(test)]

extern crate rand;
extern crate test;

use cyclotomic::fields::sparse::{random_cyclotomic, Number, is_zero};
use cyclotomic::fields::AdditiveGroup;
use cyclotomic::fields::FieldElement;
use cyclotomic::fields::MultiplicativeGroup;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use test::{black_box, Bencher};
use std::time::Instant;
use std::env;
use cyclotomic::fields::sparse::basis::convert_to_base;

fn main() {
    let gen = &mut ChaCha20Rng::seed_from_u64(12345);

    let args: Vec<String> = env::args().collect();
    let num_bench: u64 = if args.len() < 2 {
        println!("no num bench given");
        10000
    } else {
        println!("using num bench = {}", args[1]);
        args[1].parse().unwrap()
    };

    let mut nums: Vec<Number> = (1..num_bench + 1)
        .into_iter()
        .map(|_| convert_to_base(&random_cyclotomic(gen, 1000)))
        .collect();

    let start = Instant::now();
    for z in &mut nums {
        if !is_zero(&z) {
            black_box(z.mul_invert());
        }
    }

    println!("time elapsed: {} ms", start.elapsed().as_millis());
}
