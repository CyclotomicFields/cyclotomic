#![feature(test)]

extern crate antic;
extern crate rand;
extern crate test;

use antic::safe::*;
use clap::Clap;
use cyclotomic::fields::sparse::{print_gap, random_cyclotomic, Number};
use cyclotomic::fields::AdditiveGroup;
use cyclotomic::fields::MultiplicativeGroup;
use cyclotomic::fields::{Q, Z};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Result, Read};
use std::io::Write;
use std::iter::FromIterator;
use std::thread;
use std::time::Instant;
use test::black_box;
use quickcheck::Gen;

#[derive(Clap)]
#[clap(version = "1.0")]
struct Opts {
    #[clap(short, long, help = "file to output GAP benchmark code to")]
    gap_out: Option<String>,

    #[clap(short, long, default_value = "120000")]
    num_tests: usize,

    #[clap(short, long, default_value = "1")]
    threads: usize,

    #[clap(short, long, default_value = "3")]
    lower_bound_order: usize,

    #[clap(short, long, default_value = "100")]
    upper_bound_order: usize,

    #[clap(short, long, default_value = "cyclotomic", help = "cyclotomic or antic")]
    implementation: String,

    #[clap(long, default_value = "5")]
    terms: usize,

    #[clap(short, long, default_value = "100", help = "maximum absolute value of integer to use as numerator or denominator in rationals")]
    q_maximum_integer: usize,

    #[clap(long, help = "use static test data in binary")]
    static_test_data: bool
}

// Kind of like the stuff in cyclotomic::polynomial, but just for test
// data, completely unoptimized, and unsuitable for general use.
#[derive(Clone)]
struct GenericCyclotomic {
    // (exp, (numerator, denominator))
    exp_coeffs: HashMap<i64, (i64, u64)>,
    order: i64
}

// Writes a gap source file with a function f you can run that will
// do the same computations, for benchmark purposes.
fn write_gap_cycs(nums: &Vec<Number>, filename: String) -> Result<()> {
    eprintln!("writing GAP expressions to {}", filename);

    let mut file = File::create(filename).unwrap();

    file.write_all(b"SetCyclotomicsLimit(2^32-1);;\n")?;
    file.write_all(b"nop := function (x) end;;\n")?;
    file.write_all(b"f := function()\n")?;

    for i in 0..nums.len() / 6 {
        let x = &nums[6 * i].clone();
        let y = &nums[6 * i + 1].clone();
        let z = &nums[6 * i + 2].clone();
        let a = &nums[6 * i + 3].clone();
        let b = &nums[6 * i + 4].clone();
        let c = &nums[6 * i + 5].clone();

        write!(
            &mut file,
            "nop({} * {} + {} * {} + {} * {});\n",
            print_gap(x),
            print_gap(y),
            print_gap(z),
            print_gap(a),
            print_gap(b),
            print_gap(c),
        );
    }
    file.write_all(b"end;;\n")?;
    file.write_all(b"start_time := NanosecondsSinceEpoch();; f();; end_time := NanosecondsSinceEpoch();; Int((end_time-start_time)/1000000);")?;
    Ok(())
}

fn read_gap_cycs(filename: String) -> Result<Vec<GenericCyclotomic>> {
    let mut file = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    let mut cycs = vec![];

    Ok(cycs)
}

fn random_generic_cyclotomic<G>(gen: &mut G, degree: usize, num_terms: usize) -> GenericCyclotomic
where
    G: rand::RngCore,
{
    let rational_bounds = (1, 10);
    let mut result = GenericCyclotomic {
        exp_coeffs: HashMap::default(),
        order: degree as i64
    };

    for _ in 0..num_terms {
        let numerator: i64 = gen.gen_range(rational_bounds.0, rational_bounds.1) as i64;
        let denominator: u64 = gen.gen_range(rational_bounds.0, rational_bounds.1) as u64;
        let exponent: i64 = gen.gen_range(0, degree) as i64;
        result.exp_coeffs.insert(exponent, (numerator, denominator));
    }

    result
}

fn into_number(f: &GenericCyclotomic) -> Number {
    let mut coeffs = cyclotomic::fields::sparse::ExpCoeffMap::default();

    for (exp, (num, denom)) in &f.exp_coeffs {
        coeffs.insert(*exp, Q::new(Z::from(*num), Z::from(*denom)));
    }

    Number::new(f.order, &coeffs)
}

fn into_antic(
    f: &GenericCyclotomic,
    field: &mut antic::safe::NumberField,
) -> antic::safe::NumberFieldElement {
    let mut num = antic::safe::NumberFieldElement::new(field);
    for (exp, (numerator, denominator)) in &f.exp_coeffs {
        let mut term = antic::safe::NumberFieldElement::new(field);
        let mut pol = antic::safe::RationalPolynomial::new();
        let mut coeff = antic::safe::Rational::new(*numerator, *denominator);
        pol.set_coeff(*exp, &mut coeff);
        term.set_to_poly(&mut pol, field);
        let mut sum = antic::safe::NumberFieldElement::new(field);
        sum.set_to_sum_of(&mut num, &mut term, field);
        num.set(&mut sum, field);
    }
    num
}

fn main() {
    let opts: Opts = Opts::parse();
    eprintln!(
        "implementation = {}\nnum_tests = {}\nthreads = {}\nlower bound order = {}\nhigher bound order = {}",
        opts.implementation, opts.num_tests, opts.threads, opts.lower_bound_order, opts.upper_bound_order
    );

    let gen = &mut ChaCha20Rng::seed_from_u64(12345);

    // Parallelises perfectly, N threads gives N times speed.
    // We should try to beat GAP without cheating though.
    assert_eq!(opts.num_tests % opts.threads, 0);

    let num_per_test = 6;

    eprintln!("generating test data");
    let random_test_data: Vec<GenericCyclotomic> = (0..opts.num_tests * num_per_test)
        .into_iter()
        .map(|_| random_generic_cyclotomic(gen, opts.upper_bound_order, opts.terms))
        .collect();

    eprintln!("reading static test data");
    let static_test_data: Vec<GenericCyclotomic> = include!("test_data_expressions");

    let test_data = if opts.static_test_data {
        static_test_data
    } else {
        random_test_data
    };

    // cut it up into chunks for each thread
    let mut chunks: Vec<Vec<GenericCyclotomic>> = Vec::new();

    let mut start: usize = 0;
    for _ in 0..opts.threads {
        let mut chunk = vec![];
        for i in start..(start + 6 * (opts.num_tests / opts.threads)) {
            chunk.push(test_data[i].clone());
        }
        chunks.push(chunk);
        start = start + 6 * (opts.num_tests / opts.threads);
    }

    let nums: Vec<Vec<Number>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| {
            chunk
                .into_iter()
                .map(|f| into_number(&f))
                .collect()
        })
        .collect();

    let mut cyclotomic_polynomial_n =
        antic::safe::RationalPolynomial::cyclotomic(opts.upper_bound_order as u64);
    let mut cyclotomic_field_n = antic::safe::NumberField::new(&mut cyclotomic_polynomial_n);
    let antic_nums: Vec<Vec<antic::safe::NumberFieldElement>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| {
            chunk
                .into_iter()
                .map(|f| into_antic(&f, &mut cyclotomic_field_n))
                .collect()
        })
        .collect();

    if let Some(gap_out) = opts.gap_out.clone() {
        match write_gap_cycs(&nums.clone().into_iter().flatten().collect(), gap_out) {
            Ok(()) => (),
            Err(e) => eprintln!("error writing gap file: {}", e),
        }
    }

    eprintln!("starting benchmark");
    let start = Instant::now();

    if opts.implementation.as_str() == "cyclotomic" {
        number_bench(&opts, &nums)
    } else if opts.implementation.as_str() == "antic" {
        antic_bench(&opts, &mut cyclotomic_field_n, &antic_nums)
    } else {
        panic!("bad implementation {}", opts.implementation);
    }

    eprintln!("time elapsed (ms):");
    println!("{}", start.elapsed().as_millis());
}

fn number_bench(opts: &Opts, nums: &Vec<Vec<Number>>) {
    for i in 0..opts.threads {
        let num_chunk = &nums[i];
        for j in 0..opts.num_tests / opts.threads {
            let x = &mut num_chunk[6 * j].clone();
            let y = &mut num_chunk[6 * j + 1].clone();
            let z = &mut num_chunk[6 * j + 2].clone();
            let a = &mut num_chunk[6 * j + 3].clone();
            let b = &mut num_chunk[6 * j + 4].clone();
            let c = &mut num_chunk[6 * j + 5].clone();

            // black_box means don't optimise this away into a no-op
            black_box(x.mul(y).add(z.mul(a)).add(b.mul(c)));
        }
    }
}

fn antic_bench(opts: &Opts, mut cyclotomic_field_n: &mut NumberField, antic_nums: &Vec<Vec<NumberFieldElement>>) {
    for i in 0..opts.threads {
        let chunk = &antic_nums[i];
        for j in 0..opts.num_tests / opts.threads {
            let mut prod1 = antic::safe::NumberFieldElement::new(&mut cyclotomic_field_n);
            prod1.set_to_mul_of(
                &mut chunk[6 * j + 0].clone(),
                &mut chunk[6 * j + 1].clone(),
                &mut cyclotomic_field_n,
            );

            let mut prod2 = antic::safe::NumberFieldElement::new(&mut cyclotomic_field_n);
            prod2.set_to_mul_of(
                &mut chunk[6 * j + 2].clone(),
                &mut chunk[6 * j + 3].clone(),
                &mut cyclotomic_field_n,
            );

            let mut prod3 = antic::safe::NumberFieldElement::new(&mut cyclotomic_field_n);
            prod3.set_to_mul_of(
                &mut chunk[6 * j + 4].clone(),
                &mut chunk[6 * j + 5].clone(),
                &mut cyclotomic_field_n,
            );

            let mut sum1 = antic::safe::NumberFieldElement::new(&mut cyclotomic_field_n);
            sum1.set_to_sum_of(&mut prod1, &mut prod2, &mut cyclotomic_field_n);

            let mut sum2 = antic::safe::NumberFieldElement::new(&mut cyclotomic_field_n);
            sum2.set_to_sum_of(&mut sum1, &mut prod3, &mut cyclotomic_field_n);
        }
    }
}
