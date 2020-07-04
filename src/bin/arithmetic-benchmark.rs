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
use std::io::Write;
use std::iter::FromIterator;
use std::thread;
use std::time::Instant;
use test::{black_box, Bencher};
use test::NamePadding::PadNone;

fn main() {
    let gen = &mut ChaCha20Rng::seed_from_u64(12345);

    let args: Vec<String> = env::args().collect();
    let num_bench: usize = if args.len() < 2 {
        println!("no num bench given");
        120
    } else {
        println!("using num bench = {}", args[1]);
        args[1].parse().unwrap()
    };

    // Parallelises perfectly, N threads gives N times speed.
    // We should try to beat GAP without cheating, but this is our trump card.
    let num_threads: usize = 1;
    assert_eq!(num_bench % num_threads, 0);
    println!("num threads = {}", num_threads);

    let nums: Vec<Number> = (0..num_bench * 6)
        .into_iter()
        .map(|_| random_cyclotomic(gen, 100))
        .collect();

    // cut it up into chunks for each thread
    let mut chunks: Vec<Vec<Number>> = Vec::new();

    let mut start: usize = 0;
    for _ in 0..num_threads {
        let chunk = Vec::from_iter(
            nums[start..(start + 6 * (num_bench / num_threads))]
                .iter()
                .cloned(),
        );
        chunks.push(chunk);
        start = start + 6 * (num_bench / num_threads);
    }

    if args.len() > 2 {
        // writes a gap source file with a function f you can run that will
        // do the same computations, for benchmark purposes.
        let gap_out = args[2].clone();
        println!("writing to {}", gap_out);

        let mut file = File::create(gap_out).unwrap();
        file.write_all(b"f := function() return [");

        for i in 0..num_bench {
            let x = &mut nums[6 * i].clone();
            let y = &mut nums[6 * i + 1].clone();
            let z = &mut nums[6 * i + 2].clone();
            let a = &mut nums[6 * i + 3].clone();
            let b = &mut nums[6 * i + 4].clone();
            let c = &mut nums[6 * i + 5].clone();

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
        }

        file.write_all(b"]; end;");
    }

    println!("starting benchmark");
    let start = Instant::now();
    let mut threads = vec![];

    for i in 0..num_threads.clone() {
        let i = i.clone();
        let num_bench = num_bench.clone();
        let num_threads = num_threads.clone();
        let chunk = chunks[i].clone();
        let handle = thread::spawn(move || {
            for j in 0..num_bench / num_threads {
                let x = &mut chunk[6 * j].clone();
                let y = &mut chunk[6 * j + 1].clone();
                let z = &mut chunk[6 * j + 2].clone();
                let a = &mut chunk[6 * j + 3].clone();
                let b = &mut chunk[6 * j + 4].clone();
                let c = &mut chunk[6 * j + 5].clone();

                // black_box means don't optimise this away into a no-op
                black_box(x.mul(y.add(z)).eq(a.mul(b.add(c))));
            }
        });
        threads.push(handle);
    }

    for thread in threads {
        thread.join();
    }

    println!("time elapsed: {} ms", start.elapsed().as_millis());
}
