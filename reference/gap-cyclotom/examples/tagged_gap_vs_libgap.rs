use gap_cyclotom_reference::libgap::Context as LibgapContext;
use gap_cyclotom_reference::tagged::{Context as TaggedContext, Cyclotomic as TaggedCyclotomic};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn integer_terms(order: u32, density_percent: u32, salt: u32) -> Vec<(u32, i64)> {
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

fn rational_terms(order: u32, density_percent: u32, salt: u32) -> Vec<(u32, (i64, i64))> {
    let denominators = [2_i64, 3, 5, 7];
    integer_terms(order, density_percent, salt)
        .into_iter()
        .enumerate()
        .map(|(index, (exponent, numerator))| {
            (
                exponent,
                (
                    numerator,
                    denominators[(index + salt as usize) % denominators.len()],
                ),
            )
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

fn assert_matches_gap(
    gap: &LibgapContext,
    tagged: &TaggedCyclotomic,
    expected: &gap_cyclotom_reference::libgap::Cyclotomic,
) {
    let terms: Vec<_> = tagged
        .rational_terms()
        .into_iter()
        .map(|(exponent, coefficient)| {
            (
                exponent,
                (
                    coefficient.numer().to_i64().unwrap(),
                    coefficient.denom().to_i64().unwrap(),
                ),
            )
        })
        .collect();
    let actual = gap.from_terms(tagged.order(), &terms).unwrap();
    assert!(actual == *expected);
}

fn print_record(
    first: &mut bool,
    implementation: &str,
    coefficient_kind: &str,
    order: u32,
    density_percent: u32,
    terms: usize,
    iterations: u64,
    ns: f64,
) {
    if !*first {
        println!(",");
    }
    *first = false;
    print!(
        "  {{\"implementation\":\"{implementation}\",\"operation\":\"mul\",\"mode\":\"packed_gap_algorithm_tagged_coefficients\",\"coefficient_kind\":\"{coefficient_kind}\",\"order\":{order},\"density\":{:.2},\"terms_per_operand\":{terms},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
        density_percent as f64 / 100.0,
    );
}

fn main() {
    let integer_cases = [
        (5_u32, 100_u32),
        (12, 100),
        (24, 100),
        (97, 1),
        (97, 100),
        (120, 1),
        (120, 100),
        (1009, 1),
        (2520, 1),
        (10009, 1),
    ];
    let rational_cases = [
        (5_u32, 100_u32),
        (12, 100),
        (24, 100),
        (97, 1),
        (97, 100),
        (120, 1),
        (120, 100),
        (1009, 1),
    ];
    let gap = LibgapContext::new().expect("initialize unmodified libgap");
    let mut tagged = TaggedContext::new();
    let mut first = true;
    println!("[");

    for (order, density_percent) in integer_cases {
        let left_terms = integer_terms(order, density_percent, 1);
        let right_terms = integer_terms(order, density_percent, 2);
        let gap_left = gap.from_integer_terms(order, &left_terms).unwrap();
        let gap_right = gap.from_integer_terms(order, &right_terms).unwrap();
        let tagged_left = tagged.from_integer_terms(order, &left_terms).unwrap();
        let tagged_right = tagged.from_integer_terms(order, &right_terms).unwrap();
        let gap_expected = gap.mul(&gap_left, &gap_right).unwrap();
        let tagged_expected = tagged.mul(&tagged_left, &tagged_right).unwrap();
        assert_matches_gap(&gap, &tagged_expected, &gap_expected);

        let (gap_iterations, gap_ns) = measure(|| {
            black_box(gap.mul(&gap_left, &gap_right).unwrap());
        });
        let (tagged_iterations, tagged_ns) = measure(|| {
            black_box(tagged.mul(&tagged_left, &tagged_right).unwrap());
        });
        print_record(
            &mut first,
            "gap_unmodified_libgap",
            "integer",
            order,
            density_percent,
            left_terms.len(),
            gap_iterations,
            gap_ns,
        );
        print_record(
            &mut first,
            "gap_tagged_rust",
            "integer",
            order,
            density_percent,
            left_terms.len(),
            tagged_iterations,
            tagged_ns,
        );
    }

    for (order, density_percent) in rational_cases {
        let left_terms = rational_terms(order, density_percent, 1);
        let right_terms = rational_terms(order, density_percent, 2);
        let gap_left = gap.from_terms(order, &left_terms).unwrap();
        let gap_right = gap.from_terms(order, &right_terms).unwrap();
        let tagged_left = tagged.from_terms(order, &left_terms).unwrap();
        let tagged_right = tagged.from_terms(order, &right_terms).unwrap();
        let gap_expected = gap.mul(&gap_left, &gap_right).unwrap();
        let tagged_expected = tagged.mul(&tagged_left, &tagged_right).unwrap();
        assert_matches_gap(&gap, &tagged_expected, &gap_expected);

        let (gap_iterations, gap_ns) = measure(|| {
            black_box(gap.mul(&gap_left, &gap_right).unwrap());
        });
        let (tagged_iterations, tagged_ns) = measure(|| {
            black_box(tagged.mul(&tagged_left, &tagged_right).unwrap());
        });
        print_record(
            &mut first,
            "gap_unmodified_libgap",
            "rational",
            order,
            density_percent,
            left_terms.len(),
            gap_iterations,
            gap_ns,
        );
        print_record(
            &mut first,
            "gap_tagged_rust",
            "rational",
            order,
            density_percent,
            left_terms.len(),
            tagged_iterations,
            tagged_ns,
        );
    }
    println!("\n]");
}
