#![feature(test)]

extern crate antic;
extern crate rand;
extern crate test;

use antic::safe::*;
use clap::Clap;

use cyclotomic::fields::{CyclotomicFieldElement, GenericCyclotomic};

use cyclotomic::fields::structure;
use cyclotomic::fields::{dense, sparse, FieldElement};

use cyclotomic::character::inner_product;
use cyclotomic::fields::AdditiveGroupElement;
use cyclotomic::fields::MultiplicativeGroupElement;
use cyclotomic::fields::{Q, Z};

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, Result};
use std::io::{Read, Write};

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

use cyclotomic::fields::dense::basis::try_reduce;
use cyclotomic::fields::exponent::Exponent;
use cyclotomic::fields::linear_algebra::Matrix;
use cyclotomic::fields::rational::{Rational, FloatRational};
use cyclotomic::fields::simd_float;

use cyclotomic::fields::rational::FixedSizeRational;
use std::marker::PhantomData;
use num_traits::Float;

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
        about = "Read GAP cyclotomic expressions to evaluate from stdin"
    )]
    Stdin(StdinOpts),

    #[clap(
        version = "1.0",
        about = "Read some characters and do some inner products in some way"
    )]
    Character(CharacterOpts),
}

#[derive(Clap)]
struct CharacterOpts {}

#[derive(Clap)]
struct StdinOpts {
    #[clap(short, long, about = "type of each element being input")]
    element_type: String,

    #[clap(
        short,
        long,
        about = "computation to run on list of inputs",
        default_value = "sumofproducts"
    )]
    mode: String,

    #[clap(
        short,
        long,
        about = "if mode is X of Y, this is the size of chunk to run Y on"
    )]
    chunk_size: usize,
}

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
fn write_gap_cycs(nums: &Vec<sparse::Number<i64>>, filename: String) -> Result<()> {
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

fn into_structure_number(field: &structure::CyclotomicField, f: &GenericCyclotomic) -> Vec<Q> {
    let mut dense = dense::Number::from_generic(f);
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

    let big_sparse_nums: Vec<Vec<sparse::Number<Z>>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| {
            chunk
                .into_iter()
                .map(|f| sparse::Number::<Z>::from_generic(&f))
                .collect()
        })
        .collect();

    // if we're using big_order, we have to just run the benchmark for sparse,
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

    let sparse_nums: Vec<Vec<sparse::Number<i64>>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| {
            chunk
                .into_iter()
                .map(|f| sparse::Number::from_generic(&f))
                .collect()
        })
        .collect();

    let dense_nums: Vec<Vec<dense::Number>> = chunks
        .clone()
        .into_iter()
        .map(|chunk| {
            chunk
                .into_iter()
                .map(|f| dense::Number::from_generic(&f))
                .collect()
        })
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
    }

    eprintln!("starting benchmark");
    let start = Instant::now();

    // TODO: clearly should be an enum
    if top_level.implementation.as_str() == "sparse" {
        number_bench(&opts, &sparse_nums);
    } else if top_level.implementation.as_str() == "dense" {
        number_bench(&opts, &dense_nums);
    } else if top_level.implementation.as_str() == "sparse" {
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

// TODO: move all this parsing stuff into a different file

fn sexp2i64(sexp: &sexp::Sexp) -> Option<i64> {
    if let sexp::Sexp::Atom(sexp::Atom::I(order)) = sexp {
        Some(*order)
    } else {
        None
    }
}

fn sexp2z(sexp: &sexp::Sexp) -> Option<Z> {
    Some(Z::from(sexp2i64(sexp)?))
}

fn sexp2list(sexp: &sexp::Sexp) -> Option<Vec<sexp::Sexp>> {
    if let sexp::Sexp::List(sexps) = sexp {
        Some(sexps.clone())
    } else {
        None
    }
}

fn sexp2string(sexp: &sexp::Sexp) -> Option<String> {
    if let sexp::Sexp::Atom(sexp::Atom::S(string)) = sexp {
        Some(string.clone())
    } else {
        None
    }
}

macro_rules! assert_eq_return {
    ($expr1:expr, $expr2:expr) => {
        if $expr1 != $expr2 {
            return None;
        }
    };
}

fn parse_order(sexp: &sexp::Sexp) -> Option<Z> {
    let sexps = sexp2list(sexp)?;
    assert_eq_return!(sexp2string(&sexps[0])?, "order".to_owned());
    sexp2z(&sexps[1])
}

fn parse_exponent(sexp: &sexp::Sexp) -> Option<Z> {
    let sexps = sexp2list(sexp)?;
    assert_eq_return!(sexp2string(&sexps[0])?, "exponent".to_owned());
    sexp2z(&sexps[1])
}

fn parse_rational(sexp: &sexp::Sexp) -> Option<(i64, u64)> {
    let sexps = sexp2list(sexp)?;
    assert_eq_return!(sexp2string(&sexps[0])?, "rational".to_owned());
    let numerator = sexp2i64(&sexps[1])?;
    let denominator = sexp2i64(&sexps[2])?;
    Some((numerator, denominator as u64))
}

fn parse_coeffs(sexp: &sexp::Sexp) -> Option<HashMap<Z, (i64, u64)>> {
    let sexps = sexp2list(sexp)?;
    assert_eq_return!(sexp2string(&sexps[0])?, "coeffs".to_owned());

    let mut result = HashMap::new();

    for i in 1..sexps.len() {
        let coeff_sexps = sexp2list(&sexps[i])?;
        assert_eq_return!(sexp2string(&coeff_sexps[0])?, "coeff".to_string());
        let exponent = parse_exponent(&coeff_sexps[1])?;
        let (numerator, denominator) = parse_rational(&coeff_sexps[2])?;
        result.insert(exponent, (numerator, denominator));
    }

    Some(result)
}

fn sexp2cyclotomic(sexp: &sexp::Sexp) -> Option<GenericCyclotomic> {
    let sexps = sexp2list(sexp)?;
    assert_eq_return!(sexps.len(), 3);
    assert_eq_return!(sexp2string(&sexps[0])?, "cyclotomic".to_string());
    let order = parse_order(&sexps[1])?;
    let coeffs = parse_coeffs(&sexps[2])?;
    Some(GenericCyclotomic {
        order: order,
        exp_coeffs: coeffs,
    })
}

fn sexp2vector(sexp: &sexp::Sexp) -> Option<Vec<GenericCyclotomic>> {
    let mut vector = vec![];
    let cyc_list = sexp2list(sexp)?;
    assert_eq_return!(sexp2string(&cyc_list[0])?, "list".to_string());
    for i in 1..cyc_list.len() {
        vector.push(sexp2cyclotomic(&cyc_list[i])?);
    }
    Some(vector)
}

fn sexp2matrix(sexp: &sexp::Sexp) -> Option<Vec<Vec<GenericCyclotomic>>> {
    let mut matrix = vec![];
    let row_list = sexp2list(sexp)?;
    assert_eq_return!(sexp2string(&row_list[0])?, "list".to_string());
    for i in 1..row_list.len() {
        matrix.push(sexp2vector(&row_list[i])?);
    }
    Some(matrix)
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

// TODO: honestly this whole thing is horrific

fn stdin(top_level: &TopLevel, opts: &StdinOpts) {
    eprintln!("reading test matrices");
    let stdin = io::stdin();
    let lines = stdin.lock().lines();
    let sexps: Vec<String> = lines
        .into_iter()
        .map(|line| gap2sexp(line.unwrap()))
        .collect();

    if opts.element_type == "scalar".to_owned() {
        let mut scalars = vec![];
        for str in sexps {
            let sexp = sexp::parse(str.as_str()).unwrap();
            scalars.push(sexp2cyclotomic(&sexp).unwrap());
        }
        assert_eq!(scalars.len() % opts.chunk_size, 0);
        let chunks: Vec<Vec<GenericCyclotomic>> = scalars
            .chunks(opts.chunk_size)
            .map(|chunk| chunk.to_vec())
            .collect();
        let mut big_sparse_nums: Vec<Vec<sparse::Number<Z>>> = chunks
            .clone()
            .into_iter()
            .map(|chunk| {
                chunk
                    .into_iter()
                    .map(|f| sparse::Number::from_generic(&f))
                    .collect()
            })
            .collect();

        let mut sparse_nums: Vec<Vec<sparse::Number<i64>>> = chunks
            .clone()
            .into_iter()
            .map(|chunk| {
                chunk
                    .into_iter()
                    .map(|f| sparse::Number::from_generic(&f))
                    .collect()
            })
            .collect();

        let mut dense_nums: Vec<Vec<dense::Number>> = chunks
            .clone()
            .into_iter()
            .map(|chunk| {
                chunk
                    .into_iter()
                    .map(|f| dense::Number::from_generic(&f))
                    .collect()
            })
            .collect();

        eprintln!("starting scalar benchmark");
        let start = Instant::now();

        if top_level.implementation.as_str() == "sparse" {
            black_box(stdin_scalar_bench(&opts, &mut sparse_nums));
        } else if top_level.implementation.as_str() == "dense" {
            black_box(stdin_scalar_bench(&opts, &mut dense_nums));
        } else if top_level.implementation.as_str() == "sparse" {
            black_box(stdin_scalar_bench(&opts, &mut big_sparse_nums));
        } else {
            panic!("unsupported implementation!")
        }

        let elapsed = start.elapsed().as_millis();
        eprintln!("time elapsed (ms):");
        println!("{}", elapsed);
    } else if opts.element_type == "matrix".to_owned() {
        let mut matrices = vec![];
        for str in sexps {
            let sexp = sexp::parse(str.as_str()).unwrap();
            matrices.push(sexp2matrix(&sexp).unwrap());
        }
        assert_eq!(matrices.len() % opts.chunk_size, 0);
        let mut sparse_matrices: Vec<Matrix<sparse::Number<i64>, i64>> = matrices
            .clone()
            .into_iter()
            .map(|matrix| Matrix {
                value: matrix
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|generic| sparse::Number::from_generic(&generic))
                            .collect()
                    })
                    .collect(),
                exp: PhantomData,
            })
            .collect();
        eprintln!("starting matrix benchmark");
        let start = Instant::now();
        if top_level.implementation.as_str() == "sparse" {
            black_box(stdin_matrix_bench(&opts, &mut sparse_matrices));
        } else {
            panic!("bad implementation! TODO: do the others?");
        }
        let elapsed = start.elapsed().as_millis();
        eprintln!("time elapsed (ms):");
        println!("{}", elapsed);
    } else {
        panic!("bad element type {}!", opts.element_type);
    }
}

fn stdin_scalar_bench<T, E>(opts: &StdinOpts, chunks: &mut Vec<Vec<T>>)
where
    E: Exponent,
    T: CyclotomicFieldElement<E>,
{
    let mut result = T::zero_order(&E::from(1));
    for mut chunk in chunks {
        let mut chunk_result = T::one_order(&E::from(1));
        for i in 0..chunk.len() {
            chunk_result.mul(&mut chunk[i]);
        }
        result.add(&mut chunk_result);
    }
}

fn stdin_matrix_bench<T, E>(opts: &StdinOpts, matrices: &mut Vec<Matrix<T, E>>)
where
    E: Exponent,
    T: CyclotomicFieldElement<E>,
{
    let N = matrices[0].value.len();
    let mut result = Matrix::zero_matrix(N, N);

    let num_chunks = matrices.len() / opts.chunk_size;

    for i in 0..num_chunks {
        let mut chunk_result = Matrix::identity_matrix(N);
        for j in 0..opts.chunk_size {
            chunk_result = Matrix::mul(&mut chunk_result, &mut matrices[opts.chunk_size * i + j]);
        }
        result = Matrix::add(&mut result, &mut chunk_result);
    }
}

fn parse_equals(input: &str, var_name: &str) -> String {
    let toks: Vec<String> = input.split("=").map(|s| s.to_owned()).collect();
    assert_eq!(toks.len(), 2);
    assert_eq!(toks[0], var_name);
    toks[1].to_owned()
}

fn character(top_level: &TopLevel, opts: &CharacterOpts) {
    eprintln!("character product benchmark, reading from stdin");
    let mut sizes_str = String::new();
    io::stdin().read_line(&mut sizes_str);
    let sizes_nums = parse_equals(sizes_str.as_str(), "sizes")
        .split(" ")
        .skip(1)
        .map(|s| s.trim().to_owned().parse::<i64>().unwrap())
        .collect();
    //eprintln!("sizes={:?}", sizes_nums);
    let mut num_chars_str = String::new();
    io::stdin().read_line(&mut num_chars_str);
    let num_chars = parse_equals(num_chars_str.as_str(), "num_chars")
        .trim()
        .to_owned()
        .parse::<i64>()
        .unwrap();
    let mut irr_chars: Vec<Vec<GenericCyclotomic>> = vec![];
    for _ in 0..num_chars {
        let mut gap_char = String::new();
        io::stdin().read_line(&mut gap_char);
        let sexp = sexp::parse(gap2sexp(gap_char).as_str()).unwrap();
        irr_chars.push(sexp2vector(&sexp).unwrap());
    }

    // The entire rest of the input is the random character, mainly because I
    // can't work out how to get GAP to not print \ and \n everywhere
    let mut random_char_str = String::new();
    io::stdin().read_to_string(&mut random_char_str);
    let random_char_gap = parse_equals(random_char_str.as_str(), "random_char");

    let random_char_filtered = random_char_gap.replace(&['\\', '\n'][..], "");

    //eprintln!("got random_char: {}", random_char_filtered);
    let random_char: Vec<GenericCyclotomic> =
        sexp2vector(&sexp::parse(gap2sexp(random_char_filtered).as_str()).unwrap()).unwrap();

    // TODO: currently only supports sparse implementation, fix (or not)
    if top_level.implementation.as_str() == "sparse" {
        let mut sparse_irr_chars = irr_chars
            .into_iter()
            .map(|char| {
                char.into_iter()
                    .map(|f| sparse::Number::<i64>::from_generic(&f))
                    .collect()
            })
            .collect();
        let mut sparse_random_char = random_char
            .into_iter()
            .map(|f| sparse::Number::<i64>::from_generic(&f))
            .collect();
        let start = Instant::now();
        black_box(character_bench(
            &opts,
            &sizes_nums,
            &mut sparse_irr_chars,
            &mut sparse_random_char,
        ));
        let elapsed = start.elapsed().as_millis();
        eprintln!("time elapsed (ms):");
        println!("{}", elapsed);
    } else if top_level.implementation.as_str() == "simd_float" {
        let mut simd_float_irr_chars = irr_chars
            .into_iter()
            .map(|char| {
                char.into_iter()
                    .map(|f| simd_float::Number::from_generic(&f))
                    .collect()
            })
            .collect();
        let mut simd_float_random_char = random_char
            .into_iter()
            .map(|f| simd_float::Number::from_generic(&f))
            .collect();
        let start = Instant::now();
        black_box(simd_float_character_bench(
            &opts,
            &sizes_nums,
            &mut simd_float_irr_chars,
            &mut simd_float_random_char,
        ));
        let elapsed = start.elapsed().as_nanos();
        eprintln!("time elapsed (ns):");
        println!("{}", elapsed);
    } else if top_level.implementation.as_str() == "sparse_fast" {
        let mut sparse_irr_chars = irr_chars
            .into_iter()
            .map(|char| {
                char.into_iter()
                    .map(|f| sparse::Number::<i64, FixedSizeRational>::from_generic(&f))
                    .collect()
            })
            .collect();
        let mut sparse_random_char = random_char
            .into_iter()
            .map(|f| sparse::Number::<i64, FixedSizeRational>::from_generic(&f))
            .collect();
        let start = Instant::now();
        black_box(character_bench(
            &opts,
            &sizes_nums,
            &mut sparse_irr_chars,
            &mut sparse_random_char,
        ));
        let elapsed = start.elapsed().as_millis();
        eprintln!("time elapsed (ms):");
        println!("{}", elapsed);
    } else if top_level.implementation.as_str() == "sparse_float" {
        let mut sparse_irr_chars = irr_chars
            .into_iter()
            .map(|char| {
                char.into_iter()
                    .map(|f| sparse::Number::<i64, FloatRational>::from_generic(&f))
                    .collect()
            })
            .collect();
        let mut sparse_random_char = random_char
            .into_iter()
            .map(|f| sparse::Number::<i64, FloatRational>::from_generic(&f))
            .collect();
        let start = Instant::now();
        black_box(character_bench(
            &opts,
            &sizes_nums,
            &mut sparse_irr_chars,
            &mut sparse_random_char,
        ));
        let elapsed = start.elapsed().as_nanos();
        eprintln!("time elapsed (ns):");
        println!("{}", elapsed);
    } else {
        panic!("bad implementation! TODO: do the others?");
    }
}

fn simd_float_character_bench (
    opts: &CharacterOpts,
    sizes: &Vec<i64>,
    irr_chars: &mut Vec<Vec<simd_float::Number>>,
    random_char: &mut Vec<simd_float::Number>,
) {
    let mut prods: Vec<Z> = vec![];
    for irr_char in irr_chars {
        let mut prod = inner_product(sizes, irr_char, random_char);
        // We reduce because the result we're really after is the integer
        // result of the inner product - it would be cheating to not count
        // the time required for this conversion.
        prods.push(prod.nearest_int());
    }
    eprintln!("simd_float calculated prods: {:?}", prods);
}

fn character_bench<E: Exponent, Q: Rational>(
    opts: &CharacterOpts,
    sizes: &Vec<i64>,
    irr_chars: &mut Vec<Vec<sparse::Number<E, Q>>>,
    random_char: &mut Vec<sparse::Number<E, Q>>,
) {
    let mut prods: Vec<Z> = vec![];
    for irr_char in irr_chars {
        let mut prod = inner_product(sizes, irr_char, random_char);
        // We reduce because the result we're really after is the integer
        // result of the inner product - it would be cheating to not count
        // the time required for this conversion.
        prod = sparse::basis::convert_to_base(&prod);
        sparse::basis::try_reduce(&mut prod);
        prods.push(
            prod.coeffs
                .get(&E::from(0))
                .unwrap_or(&Q::zero())
                .numer()
                .clone(),
        );
    }
    eprintln!("cyclotomic calculated prods: {:?}", prods);
}

fn main() {
    let top_level: TopLevel = TopLevel::parse();

    match &top_level.subcmd {
        SubCommand::Random(random_opts) => random(&top_level, &random_opts),
        SubCommand::Stdin(stdin_opts) => stdin(&top_level, &stdin_opts),
        SubCommand::Character(character_opts) => character(&top_level, &character_opts),
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
