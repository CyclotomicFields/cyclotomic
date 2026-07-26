use cyclotomic::fields::dense;
use cyclotomic::fields::sparse::{self, ExpCoeffMap};
use cyclotomic::fields::structure::{write_dense_in_basis, CyclotomicField};
use cyclotomic::fields::MultiplicativeGroupElement;
use gap_cyclotom_reference::Context;
use std::hint::black_box;
use std::time::{Duration, Instant};

fn terms(order: u32, density_percent: u32, salt: u32) -> Vec<(u32, i64)> {
    let count = (order * density_percent).div_ceil(100).max(1);
    (0..count)
        .map(|index| {
            // 37 is coprime to every benchmark order, so this visits distinct
            // exponents before wrapping.
            let exponent = (index * 37 + salt * 11) % order;
            let magnitude = i64::from((index + salt * 3) % 5 + 1);
            let coefficient = if (index + salt) % 2 == 0 {
                magnitude
            } else {
                -magnitude
            };
            (exponent, coefficient)
        })
        .collect()
}

fn sparse_number(order: u32, terms: &[(u32, i64)]) -> sparse::Number {
    let mut coefficients = ExpCoeffMap::default();
    for &(exponent, coefficient) in terms {
        coefficients.insert(i64::from(exponent), rug::Rational::from(coefficient));
    }
    sparse::Number::new(&i64::from(order), &coefficients)
}

fn dense_number(order: u32, terms: &[(u32, i64)]) -> dense::Number {
    let mut coefficients = vec![rug::Rational::from(0); order as usize];
    for &(exponent, coefficient) in terms {
        coefficients[exponent as usize] += coefficient;
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
    // Ten cases times four implementations, with exponential calibration,
    // completes in roughly ten seconds on a development machine.
    let minimum = Duration::from_millis(100);
    let mut context = Context::new().expect("create C reference context");
    println!("[");
    let mut first = true;

    for order in [5_u32, 8, 12, 16, 24] {
        for density in [10_u32, 100] {
            let left_terms = terms(order, density, 1);
            let right_terms = terms(order, density, 2);
            let c_left = context.from_terms(order, &left_terms).unwrap();
            let c_right = context.from_terms(order, &right_terms).unwrap();
            let sparse_left = sparse_number(order, &left_terms);
            let sparse_right = sparse_number(order, &right_terms);
            let dense_left = dense_number(order, &left_terms);
            let dense_right = dense_number(order, &right_terms);
            let field = CyclotomicField::new(i64::from(order));
            let structure_left =
                write_dense_in_basis(&mut dense_number(order, &left_terms), &field.basis);
            let structure_right =
                write_dense_in_basis(&mut dense_number(order, &right_terms), &field.basis);

            let (c_iterations, c_ns) = measure(
                || {
                    black_box(context.mul(&c_left, &c_right).unwrap());
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
                ("gap_extracted_i64", c_iterations, c_ns),
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
                    "  {{\"implementation\":\"{implementation}\",\"operation\":\"mul\",\"order\":{order},\"density\":{:.2},\"terms_per_operand\":{},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
                    density as f64 / 100.0,
                    left_terms.len(),
                );
            }
        }
    }
    println!("\n]");
}
