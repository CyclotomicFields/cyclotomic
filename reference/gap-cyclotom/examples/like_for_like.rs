use cyclotomic::fields::sparse::{ExpCoeffMap, Number};
use cyclotomic::fields::MultiplicativeGroupElement;
use gap_cyclotom_reference::Context;
use std::hint::black_box;
use std::time::{Duration, Instant};

fn terms(order: u32, density_percent: u32, salt: u32) -> Vec<(u32, i64)> {
    (0..order)
        .filter(|exponent| (exponent * 37 + salt * 11 + order) % 100 < density_percent)
        .map(|exponent| (exponent, (i64::from(exponent + salt * 3) % 11) - 5))
        .filter(|term| term.1 != 0)
        .collect()
}

fn rust_number(order: u32, terms: &[(u32, i64)]) -> Number {
    let mut coefficients = ExpCoeffMap::default();
    for &(exponent, coefficient) in terms {
        coefficients.insert(i64::from(exponent), rug::Rational::from(coefficient));
    }
    Number::new(&i64::from(order), &coefficients)
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
    let minimum = Duration::from_millis(150);
    let mut context = Context::new().expect("create C reference context");
    println!("[");
    let mut first = true;

    for order in [5_u32, 8, 12, 16, 24, 30, 32, 40] {
        for density in [10_u32, 50, 100] {
            let left_terms = terms(order, density, 1);
            let right_terms = terms(order, density, 2);
            let c_left = context.from_terms(order, &left_terms).unwrap();
            let c_right = context.from_terms(order, &right_terms).unwrap();
            let rust_left = rust_number(order, &left_terms);
            let rust_right = rust_number(order, &right_terms);

            let (c_iterations, c_ns) = measure(
                || {
                    black_box(context.mul(&c_left, &c_right).unwrap());
                },
                minimum,
            );
            let (rust_iterations, rust_ns) = measure(
                || {
                    let mut left = rust_left.clone();
                    let mut right = rust_right.clone();
                    black_box(left.mul(&mut right));
                },
                minimum,
            );

            for (implementation, iterations, ns) in [
                ("gap_extracted_i64", c_iterations, c_ns),
                ("rust_sparse_rational", rust_iterations, rust_ns),
            ] {
                if !first {
                    println!(",");
                }
                first = false;
                print!(
                    "  {{\"implementation\":\"{implementation}\",\"operation\":\"mul\",\"order\":{order},\"density\":{:.2},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
                    density as f64 / 100.0
                );
            }
        }
    }
    println!("\n]");
}
