use gap_cyclotom_reference::literal::Context as RustContext;
use gap_cyclotom_reference::Context as CContext;
use std::hint::black_box;
use std::time::{Duration, Instant};

fn terms(order: u32, density_percent: u32, salt: u32) -> Vec<(u32, i64)> {
    let count = (order * density_percent).div_ceil(100).max(1);
    (0..count)
        .map(|index| {
            let exponent = (index * 37 + salt * 11) % order;
            let magnitude = i64::from((index * 13 + salt * 3) % 29 + 1);
            let coefficient = if (index + salt) % 2 == 0 {
                magnitude
            } else {
                -magnitude
            };
            (exponent, coefficient)
        })
        .collect()
}

fn measure(mut operation: impl FnMut()) -> (u64, f64) {
    let minimum = Duration::from_millis(20);
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
    let cases = [
        (5_u32, 100_u32),
        (8, 100),
        (12, 100),
        (24, 100),
        (97, 1),
        (97, 100),
        (120, 1),
        (120, 100),
        (1009, 1),
        (2520, 1),
    ];
    let mut c = CContext::new().expect("create standalone C GAP context");
    let mut rust = RustContext::new();
    let mut first = true;
    println!("[");

    for (order, density_percent) in cases {
        let left_terms = terms(order, density_percent, 1);
        let right_terms = terms(order, density_percent, 2);
        let c_left = c.from_terms(order, &left_terms).unwrap();
        let c_right = c.from_terms(order, &right_terms).unwrap();
        let rust_left = rust.from_terms(order, &left_terms).unwrap();
        let rust_right = rust.from_terms(order, &right_terms).unwrap();

        let c_expected = c.mul(&c_left, &c_right).unwrap();
        let rust_expected = rust.mul(&rust_left, &rust_right).unwrap();
        assert_eq!(rust_expected.order(), c_expected.order());
        assert_eq!(rust_expected.terms(), c_expected.terms());

        let (c_iterations, c_ns) = measure(|| {
            black_box(c.mul(&c_left, &c_right).unwrap());
        });
        let (rust_iterations, rust_ns) = measure(|| {
            black_box(rust.mul(&rust_left, &rust_right).unwrap());
        });

        for (implementation, iterations, ns) in [
            ("gap_c_translation_i64", c_iterations, c_ns),
            ("gap_literal_rust_i64", rust_iterations, rust_ns),
        ] {
            if !first {
                println!(",");
            }
            first = false;
            print!(
                "  {{\"implementation\":\"{implementation}\",\"operation\":\"mul\",\"mode\":\"literal_gap_algorithm_checked_i64\",\"order\":{order},\"density\":{:.2},\"terms_per_operand\":{},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
                density_percent as f64 / 100.0,
                left_terms.len(),
            );
        }
    }
    println!("\n]");
}
