use cyclotomic::fields::dense;
use cyclotomic::fields::sparse::{self, ExpCoeffMap};
use cyclotomic::fields::structure::{write_dense_in_basis, CyclotomicField};
use cyclotomic::fields::MultiplicativeGroupElement;
use gap_cyclotom_reference::libgap::Context;
use std::hint::black_box;
use std::time::{Duration, Instant};

type Term = (u32, (i64, i64));

#[derive(Clone, Copy)]
struct Case {
    family: &'static str,
    order: u32,
    density_percent: u32,
    coefficient_scale: i64,
    include_dense: bool,
    include_structure: bool,
}

fn terms(case: Case, salt: u32) -> Vec<Term> {
    let count = (case.order * case.density_percent).div_ceil(100).max(1);
    let denominators = [2_i64, 3, 5, 7];
    (0..count)
        .map(|index| {
            let exponent = (index * 37 + salt * 11) % case.order;
            let magnitude = case.coefficient_scale + i64::from((index + salt * 3) % 97 + 1);
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

fn measure(mut operation: impl FnMut()) -> (u64, f64) {
    let minimum = Duration::from_millis(10);
    let mut iterations = 1_u64;
    loop {
        let start = Instant::now();
        for _ in 0..iterations {
            operation();
        }
        let elapsed = start.elapsed();
        if elapsed >= minimum || iterations >= 1 << 14 {
            return (iterations, elapsed.as_nanos() as f64 / iterations as f64);
        }
        iterations *= 2;
    }
}

fn main() {
    let cases = [
        Case {
            family: "high_prime_all",
            order: 97,
            density_percent: 1,
            coefficient_scale: 1,
            include_dense: true,
            include_structure: true,
        },
        Case {
            family: "high_prime_all_large_coefficients",
            order: 97,
            density_percent: 100,
            coefficient_scale: 1_000_000_000_000,
            include_dense: true,
            include_structure: true,
        },
        Case {
            family: "highly_composite_all",
            order: 120,
            density_percent: 1,
            coefficient_scale: 1,
            include_dense: true,
            include_structure: true,
        },
        Case {
            family: "highly_composite_all_large_coefficients",
            order: 120,
            density_percent: 100,
            coefficient_scale: 1_000_000_000_000,
            include_dense: true,
            include_structure: true,
        },
        Case {
            family: "high_prime",
            order: 1009,
            density_percent: 1,
            coefficient_scale: 1,
            include_dense: true,
            include_structure: false,
        },
        Case {
            family: "high_prime_large_coefficients",
            order: 1009,
            density_percent: 100,
            coefficient_scale: 1_000_000_000_000,
            include_dense: true,
            include_structure: false,
        },
        Case {
            family: "very_divisible",
            order: 2520,
            density_percent: 1,
            coefficient_scale: 1_000_000_000_000,
            include_dense: true,
            include_structure: false,
        },
        Case {
            family: "very_divisible",
            order: 2520,
            density_percent: 100,
            coefficient_scale: 1_000_000_000_000,
            include_dense: true,
            include_structure: false,
        },
        Case {
            family: "very_high_prime_sparse",
            order: 10009,
            density_percent: 1,
            coefficient_scale: 1_000_000_000_000,
            include_dense: false,
            include_structure: false,
        },
    ];

    let context = Context::new().expect("initialize unmodified libgap");
    println!("[");
    let mut first = true;

    for case in cases {
        let left_terms = terms(case, 1);
        let right_terms = terms(case, 2);
        let gap_left = context.from_terms(case.order, &left_terms).unwrap();
        let gap_right = context.from_terms(case.order, &right_terms).unwrap();
        let sparse_left = sparse_number(case.order, &left_terms);
        let sparse_right = sparse_number(case.order, &right_terms);

        let (gap_iterations, gap_ns) = measure(|| {
            black_box(context.mul(&gap_left, &gap_right).unwrap());
        });
        let (sparse_iterations, sparse_ns) = measure(|| {
            let mut left = sparse_left.clone();
            let mut right = sparse_right.clone();
            black_box(left.mul(&mut right));
        });
        let mut measurements = vec![
            ("gap_unmodified_libgap", gap_iterations, gap_ns),
            ("rust_sparse_rational", sparse_iterations, sparse_ns),
        ];

        if case.include_dense {
            let dense_left = dense_number(case.order, &left_terms);
            let dense_right = dense_number(case.order, &right_terms);
            let (iterations, ns) = measure(|| {
                let mut left = dense_left.clone();
                let mut right = dense_right.clone();
                black_box(left.mul(&mut right));
            });
            measurements.push(("rust_dense_rational", iterations, ns));
        }

        if case.include_structure {
            let field = CyclotomicField::new(i64::from(case.order));
            let structure_left =
                write_dense_in_basis(&mut dense_number(case.order, &left_terms), &field.basis);
            let structure_right =
                write_dense_in_basis(&mut dense_number(case.order, &right_terms), &field.basis);
            let (iterations, ns) = measure(|| {
                black_box(field.mul(&structure_left, &structure_right));
            });
            measurements.push(("rust_structure_rational", iterations, ns));
        }

        for (implementation, iterations, ns) in measurements {
            if !first {
                println!(",");
            }
            first = false;
            print!(
                "  {{\"family\":\"{}\",\"implementation\":\"{implementation}\",\"operation\":\"mul\",\"coefficient_kind\":\"rational\",\"coefficient_scale\":{},\"order\":{},\"density\":{:.2},\"terms_per_operand\":{},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
                case.family,
                case.coefficient_scale,
                case.order,
                case.density_percent as f64 / 100.0,
                left_terms.len(),
            );
        }
    }
    println!("\n]");
}
