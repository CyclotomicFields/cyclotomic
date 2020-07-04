#![feature(test)]

extern crate rand;
extern crate test;

use clap::Clap;
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

#[derive(Clap)]
#[clap(version = "1.0")]
struct Opts {
    #[clap(short, long)]
    gap_out: Option<String>,

    #[clap(short, long, default_value = "120000")]
    num_tests: usize,

    #[clap(short, long, default_value = "1")]
    threads: usize,

    #[clap(short, long, default_value = "100")]
    order_max: usize,
}

fn main() {
    let opts: Opts = Opts::parse();
    eprintln!(
        "num_tests = {}\nthreads = {}\norder max = {}",
        opts.num_tests, opts.threads, opts.order_max
    );

    let gen = &mut ChaCha20Rng::seed_from_u64(12345);

    // Parallelises perfectly, N threads gives N times speed.
    // We should try to beat GAP without cheating, but this is our trump card.
    assert_eq!(opts.num_tests % opts.threads, 0);

    let nums: Vec<Number> = (0..opts.num_tests * 6)
        .into_iter()
        .map(|_| random_cyclotomic(gen, opts.order_max as i64))
        .collect();

    // cut it up into chunks for each thread
    let mut chunks: Vec<Vec<Number>> = Vec::new();

    let mut start: usize = 0;
    for _ in 0..opts.threads {
        let chunk = Vec::from_iter(
            nums[start..(start + 6 * (opts.num_tests / opts.threads))]
                .iter()
                .cloned(),
        );
        chunks.push(chunk);
        start = start + 6 * (opts.num_tests / opts.threads);
    }

    if let Some(gap_out) = opts.gap_out {
        // writes a gap source file with a function f you can run that will
        // do the same computations, for benchmark purposes.
        eprintln!("writing GAP expressions to {}", gap_out);

        let mut file = File::create(gap_out).unwrap();
        file.write_all(b"f := function() return [");

        for i in 0..opts.num_tests {
            let x = &mut nums[6 * i].clone();
            let y = &mut nums[6 * i + 1].clone();
            let z = &mut nums[6 * i + 2].clone();
            let a = &mut nums[6 * i + 3].clone();
            let b = &mut nums[6 * i + 4].clone();
            let c = &mut nums[6 * i + 5].clone();

            write!(
                &mut file,
                "{} * {} + {} * {} + {} * {},\n",
                print_gap(x),
                print_gap(y),
                print_gap(z),
                print_gap(a),
                print_gap(b),
                print_gap(c),
            );
        }

        file.write_all(b"]; end;\n");
        file.write_all(b"start_time := NanosecondsSinceEpoch();; f();; end_time := NanosecondsSinceEpoch();; Int((end_time-start_time)/1000000);");
    }

    eprintln!("starting benchmark");
    let start = Instant::now();
    let mut threads = vec![];

    for i in 0..opts.threads.clone() {
        let i = i.clone();
        let num_tests = opts.num_tests.clone();
        let num_threads = opts.threads.clone();
        let chunk = chunks[i].clone();
        let handle = thread::spawn(move || {
            for j in 0..num_tests / num_threads {
                let x = &mut chunk[6 * j].clone();
                let y = &mut chunk[6 * j + 1].clone();
                let z = &mut chunk[6 * j + 2].clone();
                let a = &mut chunk[6 * j + 3].clone();
                let b = &mut chunk[6 * j + 4].clone();
                let c = &mut chunk[6 * j + 5].clone();

                // black_box means don't optimise this away into a no-op
                black_box(x.mul(y).add(z.mul(a)).add(b.mul(c)));
            }
        });
        threads.push(handle);
    }

    for thread in threads {
        thread.join();
    }

    eprintln!("time elapsed (ms):");
    println!("{}", start.elapsed().as_millis());
}
