use cyclotomic::fields::dense;
use cyclotomic::fields::sparse::{self, ExpCoeffMap};
use cyclotomic::fields::structure::{write_dense_in_basis, CyclotomicField};
use cyclotomic::fields::MultiplicativeGroupElement;
use gap_cyclotom_reference::libgap::Context;
use std::hint::black_box;
use std::time::{Duration, Instant};

type Term = (u32, (i64, i64));

fn terms(order: u32, density_percent: u32, salt: u32) -> Vec<Term> {
    let count = (order * density_percent).div_ceil(100).max(1);
    let denominators = [2_i64, 3, 5, 7];
    (0..count)
        .map(|index| {
            let exponent = (index * 37 + salt * 11) % order;
            let magnitude = i64::from((index + salt * 3) % 5 + 1);
            let numerator = if (index + salt) % 2 == 0 {
                magnitude
            } else {
                -magnitude
            };
            let denominator = denominators[((index + salt) % 4) as usize];
            (exponent, (numerator, denominator))
        })
        .collect()
}

fn rational(coefficient: (i64, i64)) -> rug::Rational {
    rug::Rational::from((coefficient.0, coefficient.1 as u64))
}

fn sparse_number(order: u32, terms: &[Term]) -> sparse::Number {
    let mut coefficients = ExpCoeffMap::default();
    for &(exponent, coefficient) in terms {
        coefficients.insert(i64::from(exponent), rational(coefficient));
    }
    sparse::Number::new(&i64::from(order), &coefficients)
}

fn dense_number(order: u32, terms: &[Term]) -> dense::Number {
    let mut coefficients = vec![rug::Rational::from(0); order as usize];
    for &(exponent, coefficient) in terms {
        coefficients[exponent as usize] += rational(coefficient);
    }
    dense::Number::new(&i64::from(order), &coefficients)
}

fn measure(mut operation: impl FnMut(), minimum: Duration) -> (u64, f64) {
    let mut iterations = 1_u64;
    loop {
        let start = Instant::now();
        for _ in 0..iterations {
            operation();
        }
        let elapsed = start.elapsed();
        if elapsed >= minimum || iterations >= 1 << 20 {
            return (iterations, elapsed.as_nanos() as f64 / iterations as f64);
        }
        iterations *= 2;
    }
}

fn main() {
    let minimum = Duration::from_millis(100);
    let context = Context::new().expect("initialize unmodified libgap");
    println!("[");
    let mut first = true;

    for order in [5_u32, 8, 12, 16, 24] {
        for density in [10_u32, 100] {
            let left_terms = terms(order, density, 1);
            let right_terms = terms(order, density, 2);
            let gap_left = context.from_terms(order, &left_terms).unwrap();
            let gap_right = context.from_terms(order, &right_terms).unwrap();
            let sparse_left = sparse_number(order, &left_terms);
            let sparse_right = sparse_number(order, &right_terms);
            let dense_left = dense_number(order, &left_terms);
            let dense_right = dense_number(order, &right_terms);
            let field = CyclotomicField::new(i64::from(order));
            let structure_left =
                write_dense_in_basis(&mut dense_number(order, &left_terms), &field.basis);
            let structure_right =
                write_dense_in_basis(&mut dense_number(order, &right_terms), &field.basis);

            let (gap_iterations, gap_ns) = measure(
                || {
                    black_box(context.mul(&gap_left, &gap_right).unwrap());
                },
                minimum,
            );
            let (sparse_iterations, sparse_ns) = measure(
                || {
                    let mut left = sparse_left.clone();
                    let mut right = sparse_right.clone();
                    black_box(left.mul(&mut right));
                },
                minimum,
            );
            let (dense_iterations, dense_ns) = measure(
                || {
                    let mut left = dense_left.clone();
                    let mut right = dense_right.clone();
                    black_box(left.mul(&mut right));
                },
                minimum,
            );
            let (structure_iterations, structure_ns) = measure(
                || {
                    black_box(field.mul(&structure_left, &structure_right));
                },
                minimum,
            );

            for (implementation, iterations, ns) in [
                ("gap_unmodified_libgap", gap_iterations, gap_ns),
                ("rust_sparse_rational", sparse_iterations, sparse_ns),
                ("rust_dense_rational", dense_iterations, dense_ns),
                (
                    "rust_structure_rational",
                    structure_iterations,
                    structure_ns,
                ),
            ] {
                if !first {
                    println!(",");
                }
                first = false;
                print!(
                    "  {{\"implementation\":\"{implementation}\",\"operation\":\"mul\",\"coefficient_kind\":\"rational\",\"order\":{order},\"density\":{:.2},\"terms_per_operand\":{},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
                    density as f64 / 100.0,
                    left_terms.len(),
                );
            }
        }
    }
    println!("\n]");
}
