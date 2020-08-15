#![feature(test)]

extern crate antic;
extern crate rand;
extern crate test;

use antic::safe::*;
use clap::Clap;

use cyclotomic::fields::GenericCyclotomic;

use cyclotomic::fields::sparse;
use cyclotomic::fields::structure;
use cyclotomic::fields::{big_sparse, dense, FieldElement};

use cyclotomic::fields::AdditiveGroupElement;
use cyclotomic::fields::MultiplicativeGroupElement;
use cyclotomic::fields::{Q, Z};

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, Result};

use cyclotomic::fields::structure::write_dense_in_basis;
use num::Zero;
use rand::seq::IteratorRandom;
use rug::rand::RandState;
use std::cmp::min;
use std::convert::TryInto;
use std::io;
use std::process::{Command, Stdio};
use std::time::Instant;
use test::black_box;

#[derive(Clap)]
#[clap(version = "1.0")]
struct TopLevel {
    #[clap(
        short,
        long,
        default_value = "sparse",
        about = "sparse, dense, big_sparse, or antic"
    )]
    implementation: String,

    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Clap)]
enum SubCommand {
    #[clap(version = "1.0", about = "Generate random test data to benchmark")]
    Random(RandomOpts),

    #[clap(
        version = "1.0",
        about = "Read cyclotomic s-expressions to evaluate from stdin"
    )]
    Stdin(StdinOpts),
}

#[derive(Clap)]
struct StdinOpts {}

#[derive(Clap)]
struct RandomOpts {
    #[clap(short, long, about = "file to output GAP benchmark code to")]
    gap_out: Option<String>,

    #[clap(short, long, default_value = "120000")]
    num_tests: usize,

    #[clap(short, long, default_value = "1")]
    threads: usize,

    #[clap(short, long, default_value = "3")]
    lower_bound_order: usize,

    #[clap(short, long, default_value = "100")]
    upper_bound_order: usize,

    #[clap(long, default_value = "5")]
    terms: usize,

    #[clap(short, long, about = "make this fraction of the terms be nonzero")]
    density: Option<f64>,

    #[clap(short, long, about = "only works with big_sparse")]
    big_order: Option<String>,
}

// Writes a gap source file with a function f you can run that will
// do the same computations, for benchmark purposes.
fn write_gap_cycs(nums: &Vec<sparse::Number>, filename: String) -> Result<()> {
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
            sparse::print_gap(x),
            sparse::print_gap(y),
            sparse::print_gap(z),
            sparse::print_gap(a),
            sparse::print_gap(b),
            sparse::print_gap(c),
        );
    }
    file.write_all(b"end;;\n")?;
    file.write_all(b"start_time := NanosecondsSinceEpoch();; f();; end_time := NanosecondsSinceEpoch();; Int((end_time-start_time)/1000000);")?;
    Ok(())
}

fn write_gap_cycs_flat(nums: &Vec<sparse::Number>, filename: String) -> Result<()> {
    eprintln!("writing GAP expressions to {}", filename);

    let mut file = File::create(filename).unwrap();

    file.write_all(b"SetCyclotomicsLimit(2^32-1);;\n")?;
    file.write_all(b"cycs := [\n")?;

    for i in 0..nums.len() - 1 {
        let z = &nums[i];
        write!(&mut file, "{},\n", sparse::print_gap(z))?;
    }
    write!(&mut file, "{}\n", sparse::print_gap(&nums[nums.len() - 1]))?;
    file.write_all(b"];;\n")?;
    Ok(())
}

fn random_generic_cyclotomic<G>(
    gen: &mut G,
    degree: Z,
    density: Option<f64>,
    terms: usize,
) -> GenericCyclotomic
where
    G: rand::RngCore,
{
    let mut rand = RandState::new();
    let rational_bounds = (1, 10);
    let mut result = GenericCyclotomic {
        exp_coeffs: HashMap::default(),
        order: degree.clone(),
    };

    let num_terms: usize = if let Some(f) = density {
        (f * degree.to_f64()) as usize
    } else {
        terms
    };

    let mut exps = vec![];
    for _ in 0..num_terms {
        exps.push(Z::from(degree.random_below_ref(&mut rand)));
    }

    for exponent in exps {
        let numerator: i64 = gen.gen_range(rational_bounds.0, rational_bounds.1) as i64;
        let denominator: u64 = gen.gen_range(rational_bounds.0, rational_bounds.1) as u64;
        result.exp_coeffs.insert(exponent, (numerator, denominator));
    }

    result
}

fn random_generic_cyclotomics<G>(
    n: usize,
    gen: &mut G,
    degree: Z,
    density: Option<f64>,
    num_terms: usize,
) -> Vec<GenericCyclotomic>
where
    G: rand::RngCore,
{
    (0..n)
        .into_iter()
        .map(|_| random_generic_cyclotomic(gen, degree.clone(), density, num_terms))
        .collect()
}

// TODO: clearly there's some way to use the Into trait, do it

fn into_sparse_number(f: &GenericCyclotomic) -> sparse::Number {
    let mut coeffs = cyclotomic::fields::sparse::ExpCoeffMap::default();

    for (exp, (num, denom)) in &f.exp_coeffs {
        coeffs.insert(exp.to_i64().unwrap(), Q::from((*num, *denom)));
    }

    sparse::Number::new(f.order.to_i64().unwrap(), &coeffs)
}

fn into_big_sparse_number(f: &GenericCyclotomic) -> big_sparse::Number {
    let mut coeffs = cyclotomic::fields::big_sparse::ExpCoeffMap::default();

    for (exp, (num, denom)) in &f.exp_coeffs {
        coeffs.insert(exp.clone(), Q::from((*num, *denom)));
    }

    big_sparse::Number::new(&Z::from(&f.order), &coeffs)
}

fn into_dense_number(f: &GenericCyclotomic) -> dense::Number {
    let mut coeffs = vec![Q::from(0); f.order.to_usize().unwrap()];

    for (exp, (num, denom)) in &f.exp_coeffs {
        coeffs[exp.to_usize().unwrap()] = Q::from((*num, *denom));
    }

    dense::Number::new(f.order.to_i64().unwrap(), &coeffs)
}

fn into_structure_number(field: &structure::CyclotomicField, f: &GenericCyclotomic) -> Vec<Q> {
    let mut dense = into_dense_number(f);
    write_dense_in_basis(&mut dense, &field.basis)
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
        pol.set_coeff(exp.to_i64().unwrap(), &mut coeff);
        term.set_to_poly(&mut pol, field);
        let mut sum = antic::safe::NumberFieldElement::new(field);
        sum.set_to_sum_of(&mut num, &mut term, field);
        num.set(&mut sum, field);
    }
    num
}

fn random(top_level: &TopLevel, opts: &RandomOpts) {
    if opts.big_order.is_some() && top_level.implementation != "big_sparse".to_owned() {
        eprintln!("you must use big_sparse with big_order!");
        return;
    }

    eprintln!(
        "implementation = {}\nnum_tests = {}\nthreads = {}\nlower bound order = {}\nhigher bound order = {}",
        top_level.implementation, opts.num_tests, opts.threads, opts.lower_bound_order, opts.upper_bound_order
    );

    let gen = &mut ChaCha20Rng::seed_from_u64(12345);

    // Parallelises perfectly, N threads gives N times speed.
    // We should try to beat GAP without cheating though.
    assert_eq!(opts.num_tests % opts.threads, 0);

    let num_per_test = 6;

    eprintln!("generating test data");
    let test_data: Vec<GenericCyclotomic> = random_generic_cyclotomics(
        opts.num_tests * num_per_test,
        gen,
        if let Some(big_order) = &opts.big_order {
            big_order.parse::<Z>().unwrap()
        } else {
            Z::from(opts.lower_bound_order)
        },
        opts.density,
        opts.terms,
    );

    // cut it up into chunks for each thread
    let chunks: Vec<Vec<GenericCyclotomic>> = test_data
        .chunks(test_data.len() / opts.threads)
        .map(|chunk| chunk.to_vec())
        .collect();

    let big_sparse_nums: Vec<Vec<big_sparse::Number>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| {
            chunk
                .into_iter()
                .map(|f| into_big_sparse_number(&f))
                .collect()
        })
        .collect();

    // if we're using big_order, we have to just run the benchmark for big_sparse,
    // none of the implementations support big numbers
    if top_level.implementation == "big_sparse".to_owned() && opts.big_order.is_some() {
        eprintln!("running big order benchmark, only for big_sparse");
        eprintln!(
            "big order = {}",
            opts.big_order.clone().unwrap().parse::<Z>().unwrap()
        );
        let start = Instant::now();
        number_bench(&opts, &big_sparse_nums);
        eprintln!("time elapsed (ms):");
        println!("{}", start.elapsed().as_millis());
        return;
    }

    let sparse_nums: Vec<Vec<sparse::Number>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| chunk.into_iter().map(|f| into_sparse_number(&f)).collect())
        .collect();

    let dense_nums: Vec<Vec<dense::Number>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| chunk.into_iter().map(|f| into_dense_number(&f)).collect())
        .collect();

    // TODO: only supports a single fixed order, improve?
    let structure_field = structure::CyclotomicField::new(opts.lower_bound_order as i64);

    let structure_nums: Vec<Vec<Vec<Q>>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| {
            chunk
                .into_iter()
                .map(|f| into_structure_number(&structure_field, &f))
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
        match write_gap_cycs(
            &sparse_nums.clone().into_iter().flatten().collect(),
            gap_out.clone(),
        ) {
            Ok(()) => (),
            Err(e) => eprintln!("error writing gap file: {}", e),
        };
        match write_gap_cycs_flat(
            &sparse_nums.clone().into_iter().flatten().collect(),
            "flat_".to_owned() + gap_out.as_str(),
        ) {
            Ok(()) => (),
            Err(e) => eprintln!("error writing flat gap file: {}", e),
        };
    }

    eprintln!("starting benchmark");
    let start = Instant::now();

    // TODO: clearly should be an enum
    if top_level.implementation.as_str() == "sparse" {
        number_bench(&opts, &sparse_nums);
    } else if top_level.implementation.as_str() == "dense" {
        number_bench(&opts, &dense_nums);
    } else if top_level.implementation.as_str() == "big_sparse" {
        number_bench(&opts, &big_sparse_nums);
    } else if top_level.implementation.as_str() == "structure" {
        structure_number_bench(&opts, &structure_field, &structure_nums);
    } else if top_level.implementation.as_str() == "antic" {
        antic_bench(&opts, &mut cyclotomic_field_n, &antic_nums);
    } else {
        panic!("bad implementation {}", top_level.implementation);
    }

    eprintln!("time elapsed (ms):");
    println!("{}", start.elapsed().as_millis());
}

fn sexp2i64(sexp: &sexp::Sexp) -> i64 {
    if let sexp::Sexp::Atom(sexp::Atom::I(order)) = sexp {
        *order
    } else {
        panic!("couldn't parse integer")
    }
}

fn sexp2z(sexp: &sexp::Sexp) -> Z {
    Z::from(sexp2i64(sexp))
}

fn sexp2list(sexp: &sexp::Sexp) -> Vec<sexp::Sexp> {
    if let sexp::Sexp::List(sexps) = sexp {
        sexps.clone()
    } else {
        panic!("couldn't parse list")
    }
}

fn sexp2string(sexp: &sexp::Sexp) -> String {
    if let sexp::Sexp::Atom(sexp::Atom::S(string)) = sexp {
        string.clone()
    } else {
        panic!("couldn't parse string")
    }
}

fn parse_order(sexp: &sexp::Sexp) -> Z {
    let sexps = sexp2list(sexp);
    assert_eq!(sexp2string(&sexps[0]), "order".to_owned());
    sexp2z(&sexps[1])
}

fn parse_exponent(sexp: &sexp::Sexp) -> Z {
    let sexps = sexp2list(sexp);
    assert_eq!(sexp2string(&sexps[0]), "exponent".to_owned());
    sexp2z(&sexps[1])
}

fn parse_rational(sexp: &sexp::Sexp) -> (i64, u64) {
    let sexps = sexp2list(sexp);
    assert_eq!(sexp2string(&sexps[0]), "rational".to_owned());
    let numerator = sexp2i64(&sexps[1]);
    let denominator = sexp2i64(&sexps[2]);
    (numerator, denominator as u64)
}

fn parse_coeffs(sexp: &sexp::Sexp) -> HashMap<Z, (i64, u64)> {
    let sexps = sexp2list(sexp);
    assert_eq!(sexp2string(&sexps[0]), "coeffs".to_owned());

    let mut result = HashMap::new();

    for i in 1..sexps.len() {
        let coeff_sexps = sexp2list(&sexps[i]);
        assert_eq!(sexp2string(&coeff_sexps[0]), "coeff".to_string());
        let exponent = parse_exponent(&coeff_sexps[1]);
        let (numerator, denominator) = parse_rational(&coeff_sexps[2]);
        result.insert(exponent, (numerator, denominator));
    }

    result
}

fn sexp2cyclotomic(sexp: String) -> GenericCyclotomic {
    let sexp = sexp::parse(sexp.as_str()).unwrap();
    let sexps = sexp2list(&sexp);
    assert_eq!(sexps.len(), 3);
    assert_eq!(sexp2string(&sexps[0]), "cyclotomic".to_string());
    let order = parse_order(&sexps[1]);
    let coeffs = parse_coeffs(&sexps[2]);
    GenericCyclotomic {
        order: order,
        exp_coeffs: coeffs,
    }
}

fn gap2sexp(gap_cyc: String) -> String {
    let mut child = Command::new("gap2sexp")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let stdin = child.stdin.as_mut().unwrap();
    stdin.write_all(gap_cyc.as_bytes());
    let output = child.wait_with_output().unwrap();
    String::from_utf8(output.stdout).unwrap()
}

fn stdin(top_level: &TopLevel, opts: &StdinOpts) {
    let stdin = io::stdin();
    let lines = stdin.lock().lines();
    let sexps: Vec<String> = lines
        .into_iter()
        .map(|line| gap2sexp(line.unwrap()))
        .collect();
    for sexp in sexps {
        println!("{:?}", sexp2cyclotomic(sexp));
    }
}

fn main() {
    let top_level: TopLevel = TopLevel::parse();

    match &top_level.subcmd {
        SubCommand::Random(random_opts) => random(&top_level, &random_opts),
        SubCommand::Stdin(stdin_opts) => stdin(&top_level, &stdin_opts),
    }
}

fn number_bench<T>(opts: &RandomOpts, nums: &Vec<Vec<T>>)
where
    T: FieldElement + Clone,
{
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

fn structure_number_bench(
    opts: &RandomOpts,
    field: &structure::CyclotomicField,
    nums: &Vec<Vec<Vec<Q>>>,
) {
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
            black_box(field.add(
                &field.mul(x, y),
                &field.add(&field.mul(z, a), &field.mul(b, c)),
            ));
        }
    }
}

fn antic_bench(
    opts: &RandomOpts,
    mut cyclotomic_field_n: &mut NumberField,
    antic_nums: &Vec<Vec<NumberFieldElement>>,
) {
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
