#![feature(test)]

extern crate rand;
extern crate test;

use cyclotomic::fields::sparse::{print_gap, random_cyclotomic, Number};
use cyclotomic::fields::AdditiveGroup;
use cyclotomic::fields::FieldElement;
use cyclotomic::fields::MultiplicativeGroup;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::env;
use std::fmt::Debug;
use std::fs::File;
use std::time::Instant;
use test::{black_box, Bencher};
use std::io::Write;

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

    let mut nums: Vec<Number> = (1..num_bench * 6 + 1)
        .into_iter()
        .map(|_| random_cyclotomic(gen))
        .collect();

    if args.len() > 2 {
        // writes a gap source file with a function f you can run that will
        // do the same computations
        let gap_out = args[2].clone();
        println!("writing to {}", gap_out);

        let mut file = File::create(gap_out).unwrap();
        file.write_all(b"f := function() return [");

        let mut i = 0;
        for _ in 1..=num_bench {
            let x = &mut nums[i].clone();
            let y = &mut nums[i + 1].clone();
            let z = &mut nums[i + 2].clone();
            let a = &mut nums[i + 3].clone();
            let b = &mut nums[i + 4].clone();
            let c = &mut nums[i + 5].clone();

            write!(
                &mut file,
                "{} * ({} + {}) = {} * ({} + {}),\n",
                print_gap(x),
                print_gap(y),
                print_gap(z),
                print_gap(a),
                print_gap(b),
                print_gap(c),
            );

            i = i + 6;
        }

        file.write_all(b"]; end;");
    }

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
